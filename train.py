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
from mix_survival import SurvivalDataset
parser = argparse.ArgumentParser(description='Predicting survival time')

parser.add_argument('--use_cuda', '-use_cuda', default='True', type=bool, help='use cuda')
parser.add_argument('--cancer_type', '-cancer_type', default='COAD', type=str, help='cancer type choosen')
parser.add_argument('--lr', '-lr', default='5e-5', type=float, help='learning rate')
parser.add_argument('--weight_decay', '-weight_decay', default='1e-4', type=float, help='weight decay')
parser.add_argument('--momentum', '-mom', default='0.9', type=float, help='SGD momentum')
parser.add_argument('--batch_size', '-b', default='1', type=int, help='batch size')
parser.add_argument('--num_worker', '-nw', default='2', type=int, help='num_worker')
parser.add_argument('--start', '-s', default='0', type=int, help='start epoch')
parser.add_argument('--end', '-e', default='100', type=int, help='end epoch')
parser.add_argument('--experiment_id', '-eid', default='0', help='experiment id')
parser.add_argument('--experiment_name', '-name', default='prognostic_demo1', help='experiment name')
parser.add_argument('--ckpt_path_save', '-ckpt_s', default='./model_CRC_new_PLCO/', help='checkpoint path to save')
parser.add_argument('--log_path', '-lp', default='./log/', help='log path to save')
parser.add_argument('--ckpt', '-ckpt', default='./', help='checkpoint path to load')
parser.add_argument('--resume', '-r', action='store_true', help='resume from checkpoint')
parser.add_argument('--way', '-way', default='10', type=str, help='train way, 40 10 or combinate')
parser.add_argument('--load_pth_train', '-lpth_t', default='./tensor_path', help='train tensor path to load')
parser.add_argument('--load_pth_valid', '-lpth_v', default='./tensor_path', help='valid tensor path to load')
parser.add_argument('--alpha', '-a', default='1.0', type=float, help='mixup alpha')
parser.add_argument('--patch_num', '-pn', default=4, type=int, help='patch num to extract feature')
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

ckpt_path_save = os.path.join(args.ckpt_path_save, args.cancer_type + "_" + str(args.patch_num) + '_' + str(args.experiment_id))
print("--" * 20)
print(ckpt_path_save)
print("--" * 20)
if not os.path.exists(ckpt_path_save):
    os.makedirs(ckpt_path_save)
print(ckpt_path_save)

SA = SurvivalAnalysis()


# optimizer initialize
#optimizer = torch.optim.Adam(net.parameters(), lr=args.lr, betas=(0.9, 0.99), weight_decay=1e-4)
if args.optimizer == 'a':
    print('use adam')
    optimizer = torch.optim.Adam(net.parameters(), lr=args.lr, betas=(0.9, 0.99), weight_decay=args.weight_decay)


scheduler = torch.optim.lr_scheduler.StepLR(
            optimizer, step_size=2, gamma=0.9)


if args.resume:
    net = load_checkpoint(args, net)

def train(epoch, dataloader, summary):
    loss_sum = 0
    acc_sum = 0
    net.train()
    pth = ""
    length = len(dataloader)
    Prediction = torch.Tensor().cuda() #cuda(0)
    Survival = torch.Tensor().cuda()   #cuda(0)
    Observed = torch.Tensor().cuda()   #cuda(0)

    for idx, (img, T, O, _, _) in enumerate(dataloader):
        if O.sum() == 0:
            continue
        N = O.shape[0]
        img = img.cuda()
        output = net(img)

        output, T, O, at_risk, failures, ties, _ = SA.calc_at_risk(output, T, O)
        print('ties:', ties)
        T = T.cuda() #cuda(0)
        O = O.cuda() #cuda(0)
        loss = cox_cost(output, at_risk, O.reshape((N, 1)), failures, ties)
        loss.register_hook(lambda g: print(g))

        optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(net.parameters(), 10)
        optimizer.step()

        Prediction = torch.cat((Prediction, output))
        Survival = torch.cat((Survival, T.float()))
        Observed = torch.cat((Observed, O.float()))
    scheduler.step()
    Prediction, Survival, Observed, at_risk, failures, ties, _ = SA.calc_at_risk(Prediction, Survival.cpu(), Observed.cpu())

    CI = concordance_index(Survival.cpu().detach().numpy(), -Prediction.cpu().detach().numpy(),
                           Observed.cpu().detach().numpy())
    loss = cox_cost(Prediction, at_risk, Observed.reshape((Observed.shape[0],1)).cuda(), failures, ties)
    print("loss:", loss.item(), "CI:", CI.item())
    summary['loss'] = loss.item()
    summary['CI'] = CI.item()
    summary['lr'] = optimizer.param_groups[0]['lr']
    return summary

def valid(dataloader, summary):
    net.eval()
    length = len(dataloader)
    Prediction = torch.Tensor().cuda()
    Survival = torch.Tensor().cuda()
    Observed = torch.Tensor().cuda()

    with torch.no_grad():
        for idx, (img, T, O, _, _ ) in enumerate(dataloader):
            img = img.cuda()
            N = O.shape[0]
            img = img.cuda()
            output = net(img)

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
    loss = cox_cost(Prediction, at_risk, Observed.reshape((Observed.shape[0],1)).cuda(), failures, ties)
    print("loss:", loss.item(), "CI:", CI.item())
    summary['loss'] = loss.item()
    summary['CI'] = CI.item()
    return summary


if __name__ == '__main__':
    pre = str(args.experiment_id)
    patch_num = args.patch_num
    
    train_file_paths = ['./datasets/files/1030_PLCO.tsv']
    train_data_paths = ['../TLS_Quantitative_Analysis/TLS_Segmantation/predictions/PLCO_Predictions/']
    train_tumor_paths = ['../datasets/PLCO']
    
    train_data = SurvivalDataset(train_file_paths, train_data_paths, train_tumor_paths, way="train", val=False, tumor_ratio=0.001, patch_num=patch_num)
    # val_data = SurvivalDataset(train_file_paths, train_data_paths, train_tumor_paths, way="train", val=False, tumor_ratio=0.001, patch_num=patch_num)
    #val_data = SurvivalDataset(train_file_paths, train_data_paths, train_tumor_paths, way="val", val=False)

    batch_size_train = 4
    num_workers = 2
    from tqdm import tqdm
    train_dataloader = torch.utils.data.DataLoader(train_data,
                                               batch_size=batch_size_train,
                                               num_workers=num_workers,
                                               drop_last=True,
                                               shuffle=True)

    batch_size_val = 4
    num_workers = 2

    val_dataloader = torch.utils.data.DataLoader(val_data,
                                               batch_size=batch_size_val,
                                               num_workers=num_workers,
                                               drop_last=False,
                                               shuffle=False)



    summary_train = {'epoch': 0, 'fp': 0, 'tp': 0, 'Neg': 0, 'Pos': 0}
    summary_valid = {'loss': float('inf'), 'acc': 0}
    summary_writer = SummaryWriter(log_path)
    loss_valid_best = float('inf')
    for epoch in range(args.start, args.end):

        summary_train = train(epoch, train_dataloader, summary_train)

        summary_writer.add_scalar(
            'train/loss', summary_train['loss'], epoch)
        summary_writer.add_scalar(
            'train/CI', summary_train['CI'], epoch)
        if epoch % 1 == 0:
            torch.save({'epoch': summary_train['epoch'],
                        'state_dict': net.state_dict()},
                       (ckpt_path_save + '/' + str(epoch) + '.ckpt'))

        summary_valid = valid(val_dataloader, summary_valid)

        summary_writer.add_scalar(
            'valid/loss', summary_valid['loss'], epoch)
        summary_writer.add_scalar(
            'valid/CI', summary_valid['CI'], epoch)
        summary_writer.add_scalar(
            'learning_rate', summary_train['lr'], epoch
        )
        print('train/loss', summary_train['loss'], epoch)
        print('train/CI', summary_train['CI'], epoch)
        print('valid/loss', float(summary_valid['loss']), epoch)
        print('valid/CI', summary_valid['CI'], epoch)

        if summary_valid['loss'] < loss_valid_best:
            loss_valid_best = summary_valid['loss']
            print("loss_valid_best is", loss_valid_best)
            torch.save({'epoch': summary_train['epoch'],
                        'optimizer': optimizer.state_dict(),
                        'state_dict': net.state_dict()},
                       os.path.join(ckpt_path_save, 'best.ckpt'))
    summary_writer.close()
