import argparse
import os
os.environ['CUDA_VISIBLE_DEVICES']='1'
import sys
import torch
import torch.utils.data
import torch.backends.cudnn as cudnn
from tensorboardX import SummaryWriter
from utils.Survival_Aanlysis import SurvivalAnalysis
from utils.RiskLayer import cox_cost
import random
from Prognostic import RiskModel
from lifelines.utils import concordance_index
import numpy as np
#from survival import SurvivalDataset
from mix_survival import SurvivalDataset
import glob
import pandas as pd
parser = argparse.ArgumentParser(description='Predicting survival time')

parser.add_argument('--use_cuda', '-use_cuda', default='True', type=bool, help='use cuda')
parser.add_argument('--cancer_type', '-cancer_type', default='COAD', type=str, help='cancer type choosen')
parser.add_argument('--lr', '-lr', default='5e-5', type=float, help='learning rate')
parser.add_argument('--weight_decay', '-weight_decay', default='1e-4', type=float, help='weight decay')
parser.add_argument('--momentum', '-mom', default='0.9', type=float, help='SGD momentum')
parser.add_argument('--batch_size', '-b', default='1', type=int, help='batch size')
parser.add_argument('--num_patch', '-n', default=24, type=int, help='number patches')
parser.add_argument('--num_worker', '-nw', default='2', type=int, help='num_worker')
parser.add_argument('--start', '-s', default='0', type=int, help='start epoch')
parser.add_argument('--end', '-e', default='100', type=int, help='end epoch')
parser.add_argument('--experiment_id', '-eid', default='0', help='experiment id')
parser.add_argument('--experiment_name', '-name', default='prognostic_demo1', help='experiment name')
parser.add_argument('--ckpt_path_save', '-ckpt_s', default='./model/', help='checkpoint path to save')
parser.add_argument('--log_path', '-lp', default='./log/', help='log path to save')
parser.add_argument('--ckpt', '-ckpt', default='model/COAD_5/best.ckpt', help='checkpoint path to load')
parser.add_argument('--resume', '-r', action='store_true', help='resume from checkpoint')
parser.add_argument('--way', '-way', default='10', type=str, help='train way, 40 10 or combinate')
parser.add_argument('--load_pth_train', '-lpth_t', default='./tensor_path', help='train tensor path to load')
parser.add_argument('--load_pth_valid', '-lpth_v', default='./tensor_path', help='valid tensor path to load')
parser.add_argument('--alpha', '-a', default='1.0', type=float, help='mixup alpha')
parser.add_argument('--device_ids', default='1', type=str, help='comma separated indices of GPU to use,'
                                                                      ' e.g. 0,1 for using GPU_0'
                                                                      ' and GPU_1, default 0.')
parser.add_argument('--drop_group', '-drop_group', default='3,4', help='drop groups')
parser.add_argument('--drop_prob', '-drop_prob', default='0.1', type=float, help='drop prob')
parser.add_argument('--freeze', '-f', action='store_true', help='Freeze convolutional layer parameters')
parser.add_argument('--type-key', '-type-key', default='tumor', type=str, help='tumor or tumor_beside or fibrous_tissue')
parser.add_argument('--experimentway', '-eway', default='prognosis', type=str, help='prognosis or replase')
parser.add_argument('--use_std', '-std', default='use', type=str, help='use std as feature, u:use, o:only, n:not use ')
parser.add_argument('--optimizer', '-o', default='a', type=str, help='choose optimizer:a(adam), s(sgd), '
                                                                     'Adadelta(Adadelta), m(MinimalLaycaSGD) '
                                                                     'or l(LaycaSGD)')
args = parser.parse_args()
cudnn.benchmark = True

def seed_torch(seed=42):
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True
seed_torch()



def load_checkpoint(args, net):
    print("Use ckpt: ", args.ckpt)
    assert len(args.ckpt) != 0, "Please input a valid ckpt_path"
    checkpoint = torch.load(args.ckpt)
    pretrained_dict = checkpoint['state_dict']
    net.load_state_dict(pretrained_dict)
    return net


# initialize model
net = RiskModel().cuda()
net = torch.nn.DataParallel(net)



log_path = os.path.join(args.log_path, args.experiment_name + "_" + str(args.experiment_id))
if not os.path.isdir(log_path):
    os.makedirs(log_path)

SA = SurvivalAnalysis()



def valid(dataloader, summary):
    infos = {}
    net.eval()
    #ref_dir = '/home/users/zqchen/datasets/TCGA_{}/'.format(args.cancer_type)
    if args.cancer_type != 'PLCO':
        ref_dir = '/home/users/zqchen/datasets/TCGA_{}/'.format(args.cancer_type)
    else:
        ref_dir = '/home/users/zqchen/datasets/{}/'.format(args.cancer_type)
    length = len(dataloader)
    Prediction = torch.Tensor().cuda()
    Survival = torch.Tensor().cuda()
    Observed = torch.Tensor().cuda()

    with torch.no_grad():
        for idx, (img, T, O, _, names ) in enumerate(dataloader):
            img = img.cuda()
            N = O.shape[0]
            img = img.cuda()
            output = net(img)
            for name, score in zip(names, output):
                f = glob.glob(os.path.join(ref_dir, name.split('/')[-1], '*svs'))[0]
                if 'PLCO' in f:
                    name = f.split('/')[-2]
                else:
                    names = f.split('/')[-1].split('-')
                    name = names[0] + '-' + names[1] + '-' + names[2]
                
                infos[name] = score.cpu().numpy()[0]

            output, T, O, at_risk, failures, ties, _ = SA.calc_at_risk(output, T, O)

            T = T.cuda()
            O = O.cuda()

            loss = cox_cost(output, at_risk, O.reshape((N, 1)), failures, ties)
            #print("loss:", loss.item())
            Prediction = torch.cat((Prediction, output))
            Survival = torch.cat((Survival, T.float()))
            Observed = torch.cat((Observed, O.float()))

    Prediction, Survival, Observed, at_risk, failures, ties, _ = SA.calc_at_risk(Prediction, Survival.cpu(), Observed.cpu())

    
    
    CI = concordance_index(Survival.cpu().detach().numpy(), -Prediction.cpu().detach().numpy(),
                           Observed.cpu().detach().numpy())
    summary['CI'] = CI.item()
    return summary, infos


if __name__ == '__main__':
    pre = str(args.experiment_id)
    if args.cancer_type != 'PLCO':
        val_file_path = './datasets/files/filter_TCGA-{}.survival.tsv'.format(args.cancer_type)
        tumor_path = '../datasets/TCGA_{}/'.format(args.cancer_type)
    else:
        val_file_path = './datasets/files/{}.survival.tsv'.format(args.cancer_type)
        tumor_path = '../datasets/{}/'.format(args.cancer_type)
    data_path = './Tumor_Segmantation/predictions/{}_Predictions/'.format(args.cancer_type)
    
    
    
    if args.resume:
        print('Load model ------>', args.ckpt)
        
        net = load_checkpoint(args, net)

    #val_data = SurvivalDataset(val_file_path, data_path, tumor_path, way="val", val=True, patch_num=int(args.num_patch))
    val_data = SurvivalDataset([val_file_path], [data_path], [tumor_path], way="train", val=False, patch_num=int(args.num_patch))

    from tqdm import tqdm

    batch_size_val = 4
    num_workers = 2

    val_dataloader = torch.utils.data.DataLoader(val_data,
                                               batch_size=batch_size_val,
                                               num_workers=num_workers,
                                               drop_last=False,
                                               shuffle=False)
    data = pd.read_table(val_file_path)
    Survival = data.time.values
    Observed = data.OS.values
    Id =list(data._PATIENT.values)
    
    #ckpts = os.listdir('model/PLCO_mix_check_32_1018')
    
    #exists = list(str(i) + '.ckpt' for i in range(1, 31))
    
    #for exist in exists:
    #    ckpts.remove(exist)

    #for ckpt in ckpts:
    #    checkpoint = torch.load(os.path.join('model/PLCO_mix_check_32_1018', ckpt))
    #    pretrained_dict = checkpoint['state_dict']
    #    print("load ckpt is ==>", ckpt)
    #    net.load_state_dict(pretrained_dict)

    summary = {} 
    summary, infos = valid(val_dataloader, summary) 
    
    
    names = []
    survivals= []
    observeds=[]
    risk_scores = []
    
    for name, survival, observed in zip(Id, Survival, Observed):
        if '-11A' in name:
            continue
        if name in infos:
            names.append(name)
            survivals.append(survival)
            observeds.append(observed)
            risk_scores.append(infos[name])

    print(args.ckpt, summary)
    #dic = {'name': names, 'time': survivals, 'OS': observeds, 'risk_scores': risk_scores}
    #df = pd.DataFrame(dic)
    #df.to_csv('confirm_csv/{}_{}_risk.csv'.format(args.num_patch, args.cancer_type))

