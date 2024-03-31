import torch
import torch.nn as nn
import numpy as np
import os
from PIL import Image
from PIL import ImageFile
import xlrd
import pandas as pd
import openslide
ImageFile.LOAD_TRUNCATED_IMAGES = True
import numpy as np
from torchvision import transforms
import random
import cv2
import glob
random.seed(0)

# df.OS.values
class SurvivalDataset(nn.Module):
    def __init__(self, file_path, data_path, tumor_path, normalize=True, way="train", val=False, tumor_ratio=0.5, patch_num=2):
        super(SurvivalDataset, self).__init__()
        self.transform = transforms.Compose([
            transforms.ToTensor(),
            transforms.ToPILImage()
        ])
        self._color_jitter = transforms.ColorJitter(64.0 / 255, 0.75, 0.25, 0.04)
        self.tumor_path = tumor_path
        self._normalize = normalize
        self._crop_size = 1024
        self._way = way
        self.val = val
        self.data_path = data_path
        self.Observed = []
        self.Survival = []
        self.file_path = file_path
        self.img_path = []
        self.annos = self.read_table()
        self.keys = self.process_dirs()
        self.tumor_ratio = tumor_ratio
        self.patch_num = patch_num
        
        
        print(self.keys, self.tumor_path)
        
        for k in self.keys:
            file_path = glob.glob(os.path.join(self.tumor_path, k, "*svs"))[0]
            keys = file_path.split('/')[-1].split('-')
            key = keys[0] + '-' + keys[1] + '-' + keys[2]
            if key not in self.annos.keys():
                continue
            
            with open(os.path.join(self.data_path, k, 'tumor.txt'), 'r') as f:
                lines = f.readlines()
                length = 0
                for line in lines:
                    content = line.strip().split('\t')
                    ratio = float(content[1])
                    if ratio > self.tumor_ratio:
                        length += 1
                #length = len(os.listdir(os.path.join(self.data_path, k, 'intra-tumors')))
                if length >= self.patch_num:
                    self.img_path.append(os.path.join(self.data_path, k))
                    self.Observed.append(self.annos[key][1])
                    self.Survival.append(self.annos[key][0])
        

    def read_table(self):
        data = pd.read_table(self.file_path)
        Survival = data.time.values
        Observed = data.OS.values
        Id =list(data._PATIENT.values)
        return dict(zip(Id, zip(Survival, Observed)))

    def process_dirs(self):
        dirs = os.listdir(self.data_path)
        keys_list = []
        for d in dirs:
            try:
                file_path = glob.glob(os.path.join(self.tumor_path, d, "*svs"))[0]
                keys = file_path.split('/')[-1].split('-')
                key = keys[0] + '-' + keys[1] + '-' + keys[2]
            except:
                continue
            if key in self.annos.keys():
                keys_list.append(d)
        return keys_list
    
    def __len__(self):
        return len(self.img_path)
    
    def _preprocess(self, im):
        #if self._way == "train" or self._way == "val":
        if self._way == "train":
            #im = self._color_jitter(im)
            if np.random.rand() > 0.5:
                im = im.transpose(Image.FLIP_LEFT_RIGHT)
            num_rotate = np.random.randint(0, 4)
            im = im.rotate(90 * num_rotate)
            # crop
            im = np.array(im)
            img_h, img_w, _ = im.shape
            pad_h = max(self._crop_size - img_h, 0)
            pad_w = max(self._crop_size - img_w, 0)
            if pad_h > 0 or pad_w > 0:
                img_pad = cv2.copyMakeBorder(im, 0, pad_h, 0,
                                             pad_w, cv2.BORDER_CONSTANT,
                                             value=(255.0, 255.0, 255.0))
            else:
                img_pad = im
            h_off = random.randint(0, img_h - self._crop_size)
            w_off = random.randint(0, img_w - self._crop_size)
            im = np.asarray(img_pad[h_off: h_off + self._crop_size, w_off: w_off + self._crop_size])
        im = np.array(im, dtype=np.float32).transpose((2, 0, 1))
        if self._normalize:
            im = (im - 128.0) / 128.0
        return im

    def __getitem__(self, idx):
        if self._way == 'train':
            ExtractSize = int(self._crop_size * 1.5)
            #ExtractSize = self._crop_size
        elif self._way == 'val':
            ExtractSize = self._crop_size
            #ExtractSize = int(self._crop_size * 1.5)
        
        RepetitionRate = 0
        T = self.Survival[idx]
        O = self.Observed[idx]
        d = self.img_path[idx].split('/')[-1]
        # ../immune_therapy/algorithms/TLS_Segmantation/predictions/TCGA_STAD_Predictions/08db7ebc-bbfa-4d7e-9e63-b88100d6693e/intra-tumors/
        #files = os.listdir(self.img_path[idx] + '/' + 'intra-tumors')
        files = []
        ratios = []
        with open(os.path.join(self.img_path[idx], 'tumor.txt'), 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip()
                contents = line.split('\t')
                ratio = float(contents[1])
                name = contents[0]
                if ratio > self.tumor_ratio:
                    files.append(name)
                    ratios.append(ratio)

        zipped = zip(ratios, files)
        sorted_queue1 = sorted(zipped, key=lambda x: x[0])
        files = [item[1] for item in sorted_queue1]


        curr_d = os.path.join(self.tumor_path, d)
        ref_img_path = glob.glob(os.path.join(curr_d, "*svs"))[0]
        slide = openslide.OpenSlide(ref_img_path)
        pth = []
        if self.val:
            sample_num = self.patch_num
        else:
            sample_num = self.patch_num
        if len(files) < sample_num:
            imgs = files
        else:
            imgs = files[len(files)-sample_num:]#random.sample(files, sample_num)
        choose_name = self.img_path[idx]
        #count = len(imgs)
        Imgs = []
        for img in imgs:
            img_name = img.replace('.png', '')
            j_i = img_name.split('_')
            j = int(float(j_i[0]))
            i = int(float(j_i[1]))
            Img = slide.read_region((int(j * ExtractSize * (1 - RepetitionRate)), int(i * ExtractSize * (1 - RepetitionRate))), 0, (ExtractSize, ExtractSize)).convert('RGB')
            pth.append(img)

            #Img = Image.open(self.img_path[idx] + '/images/' + img)
            Img1 = self._preprocess(Img)
            Imgs.append(Img1)
            #Img2 = self._preprocess(Img)
            #Imgs.append(Img2)

        Imgs = np.array(Imgs)
        return Imgs, T, O, pth, choose_name
