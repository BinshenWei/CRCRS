import numpy as np
from torch import nn
import torch
from torchvision import models, transforms, datasets
import torch.nn.functional as F
import pretrainedmodels
from nystrom_attention import NystromAttention
from utils.distribution_pooling_filter import DistributionPoolingFilter

class FeatureExtractor(nn.Module):
    def __init__(self, num_in=1024, num_out=32):
        super(FeatureExtractor, self).__init__()
        self.num_ftrs = num_in
        self.num_out = num_out
        self._model_conv.fc = nn.Linear(self.num_ftrs, self.num_out)
        self.sigmoid = nn.Sigmoid()
    def forward(self, x):
        out = self._model_conv(x)
        out = self.sigmoid(out)
        return out

class RepresentationTransformation(nn.Module):
    def __init__(self, num_in=32, num_out=1):
        super(RepresentationTransformation, self).__init__()
        self.fc = nn.Sequential(
                nn.Dropout(0.5),
                nn.Linear(num_in, 128),
                nn.ReLU(),
                nn.Dropout(0.5),
                nn.Linear(128, 32),
                nn.ReLU(),
                nn.Dropout(0.5),
                nn.Linear(32, num_out)
                    )

    def forward(self, x):
        out = self.fc(x)
        return out

class Attention(nn.Module):
    def __init__(self, num_in=32, num_instances=32):
        super(Attention, self).__init__()
        self._num_instances = num_instances
        self.fc = nn.Sequential(
                        nn.Linear(num_in, 128),
                        nn.Tanh(),
                        nn.Linear(128, 1)
                        )

    def forward(self, x):
        out = self.fc(x)
        out = torch.reshape(out,(-1,self._num_instances,1))
        out = F.softmax(out, dim=1)
        return out

class Attention2(nn.Module):
    def __init__(self, num_in=32, num_instances=32):
        super(Attention2, self).__init__()
        self._num_instances = num_instances
        self.fc = nn.Sequential(
                        nn.Linear(num_in, 128),
                        nn.ReLU(),
                        nn.Linear(128, 1)
                    )

    def forward(self, x):
        out = self.fc(x)
        out = torch.reshape(out,(-1,self._num_instances,1))
        out = torch.sigmoid(out)
        return out

class TransLayer(nn.Module):

    def __init__(self, norm_layer=nn.LayerNorm, dim=512):
        super().__init__()
        self.norm = norm_layer(dim)
        self.attn = NystromAttention(
            dim = dim,
            dim_head = dim//8,
            heads = 8,
            num_landmarks = dim//2,    # number of landmarks
            pinv_iterations = 6,    # number of moore-penrose iterations for approximating pinverse. 6 was recommended by the paper
            residual = True,         # whether to do an extra residual with the value or not. supposedly faster convergence if turned on
            dropout=0.1
        )

    def forward(self, x):
        x = x + self.attn(self.norm(x))
        return x

class PPEG(nn.Module):
    def __init__(self, dim=512):
        super(PPEG, self).__init__()
        self.proj = nn.Conv2d(dim, dim, 7, 1, 7//2, groups=dim)
        self.proj1 = nn.Conv2d(dim, dim, 5, 1, 5//2, groups=dim)
        self.proj2 = nn.Conv2d(dim, dim, 3, 1, 3//2, groups=dim)

    def forward(self, x, H, W):
        B, _, C = x.shape
        cls_token, feat_token = x[:, 0], x[:, 1:]
        cnn_feat = feat_token.transpose(1, 2).view(B, C, H, W)
        x = self.proj(cnn_feat)+cnn_feat+self.proj1(cnn_feat)+self.proj2(cnn_feat)
        x = x.flatten(2).transpose(1, 2)
        x = torch.cat((cls_token.unsqueeze(1), x), dim=1)
        return x




class RiskModel(nn.Module):
    def __init__(self, backbone='resnet50', num_classes=1):
        super(RiskModel, self).__init__()
        self.num_classes = num_classes
        #self.way = way
        self.backbone_arch = backbone
        if self.backbone_arch == 'resnet50':
            self.model = models.resnet50()
            self.model.load_state_dict(torch.load("/home/users/zqchen/.cache/torch/hub/checkpoints/resnet50-19c8e357.pth"))

        self.avgpool = nn.AdaptiveAvgPool2d(output_size=1)
        self.fc = nn.Linear(1024, self.num_classes, bias=False)

        self.conv1 = nn.Conv2d(3,64,3,padding=1)
        self.bn1 = nn.BatchNorm2d(64)
        self.relu = nn.ReLU(inplace=True)
        
        #stage 1
        self.encoder1 = self.model.layer1
        #stage 2
        self.encoder2 = self.model.layer2
        #stage 3
        self.encoder3 = self.model.layer3
        #stage 4
        #self.encoder4 = self.model.layer4
        self.avgpool = nn.AdaptiveAvgPool2d(output_size=1)
        self.maxpool = nn.MaxPool2d(kernel_size=3, stride=2, padding=1)
        self.num_features = 1024
        self.num_instances = 2
        self.num_bins=11
        self.sigma=0.1
        self._attention = Attention(num_in=self.num_features, num_instances=self.num_instances)
        self._attention2 = Attention2(num_in=self.num_features, num_instances=self.num_instances)
        self._distribution_pooling_filter = DistributionPoolingFilter(num_bins=self.num_bins, sigma=self.sigma)
        self._representation_transformation = RepresentationTransformation(num_in=self.num_features*self.num_bins, num_out=self.num_classes)
        



        self.pos_layer = PPEG(dim=512)
        self._fc1 = nn.Sequential(nn.Linear(1024, 512), nn.ReLU())
        self.cls_token = nn.Parameter(torch.randn(1, 1, 512))
        self.layer1 = TransLayer(dim=512)
        self.layer2 = TransLayer(dim=512)
        self.norm = nn.LayerNorm(512)
        self._fc2 = nn.Linear(512, self.num_classes)

    
    def forward(self, x):
        batch_size,num_patch,_,img_size,_ = x.shape
        x = x.view(-1,3,img_size,img_size)
        
        x = self.conv1(x)
        x = self.bn1(x)
        x = self.relu(x)
        x = self.maxpool(x)

        x = self.encoder1(x)
        x = self.encoder2(x)
        x = self.encoder3(x)
        
        #x = self.encoder4(x)
        logits = x.view(batch_size, num_patch, x.size(1), x.size(-1), x.size(-1))
        '''
        extracted_features = torch.squeeze(torch.squeeze(self.avgpool(logits), dim=-1), dim=-1)
        attention_values = self._attention(extracted_features)
        attention_values2 = self._attention2(extracted_features)
        extracted_features = torch.reshape(extracted_features,(-1,self.num_instances,self.num_features))
        out = attention_values2*extracted_features
        out = self._distribution_pooling_filter(out, attention_values)
        out = torch.reshape(out,(-1, self.num_features*self.num_bins))
        out = self._representation_transformation(out)
        return out
        x, _ = logits.max(1)
        x = x.squeeze(1)
        x = self.avgpool(x)
        x = x.view(x.size(0), -1)
        x = self.fc(x)

        return x
        '''
        
        
        #logits = x.view(batch_size, num_patch, x.size(1), x.size(-1), x.size(-1))
        h = torch.squeeze(torch.squeeze(self.avgpool(logits), dim=-1), dim=-1)
        h = self._fc1(h) #[B, n, 512]
        #---->pad
        H = h.shape[1]
        _H, _W = int(np.ceil(np.sqrt(H))), int(np.ceil(np.sqrt(H)))
        add_length = _H * _W - H
        h = torch.cat([h, h[:,:add_length,:]],dim = 1) #[B, N, 512]

        #---->cls_token
        B = h.shape[0]
        cls_tokens = self.cls_token.expand(B, -1, -1).cuda()
        h = torch.cat((cls_tokens, h), dim=1)

        #---->Translayer x1
        h = self.layer1(h) #[B, N, 512]

        #---->PPEG
        h = self.pos_layer(h, _H, _W) #[B, N, 512]

        #---->Translayer x2
        h = self.layer2(h) #[B, N, 512]

        #---->cls_token
        h = self.norm(h)[:,0]
        #---->predict
        x = self._fc2(h) #[B, n_classes]
        
        print("logits shape is ", logits.shape)
        return x, logits
        
