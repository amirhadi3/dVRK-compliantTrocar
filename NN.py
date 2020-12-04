import torch
from torch.autograd import Variable
import torch.nn.functional as F
import torch.utils.data as Data

import matplotlib.pyplot as plt
#matplotlib inline

import numpy as np
import imageio
from numpy import genfromtxt


# settings
TRAIN = True

# read the training data
repo = 'ingripper calibration/fast - 5 cycles - pi_3/'
#repo = 'ingripper calibration/fast - 5 cycles - pi_3/'
ati = np.genfromtxt(repo + 'atiData.txt', delimiter=' ',dtype=float)
diffsig = np.genfromtxt(repo + 'diffData.txt', delimiter=' ',dtype=float)
sumsig = np.genfromtxt(repo + 'sumData.txt', delimiter=' ',dtype=float)
pos = np.genfromtxt(repo + 'posData.txt', delimiter=' ',dtype=float)
vel = np.genfromtxt(repo + 'velData.txt', delimiter=' ',dtype=float)
eff = np.genfromtxt(repo + 'effData.txt', delimiter=' ',dtype=float)
#temp = np.genfromtxt(repo + 'tempData.txt', delimiter=' ',dtype=float)
#temp = temp[:,np.newaxis]

# build the feature vector
diffsig = diffsig-np.mean(diffsig[0:500,:],axis = 0)
ati = ati-np.mean(ati[0:500,:],axis = 0)
Nsig = diffsig/sumsig
Nsig = np.hstack((Nsig,Nsig**2))

x = np.hstack((Nsig,pos,vel,eff));
mean_x = np.mean(x,axis=0)
stddev_x = np.std(x,axis=0)
x = (x-mean_x)/stddev_x
max_ati = np.max(np.abs(ati),axis=0);
y = ati/max_ati

xss = x[0:x.shape[0]:50,:]
yss = y[0:y.shape[0]:50,:]

# split data to training and test sets
dataSize = xss.shape[0]
trainSize = int(np.floor(0.7*dataSize))
testSize = dataSize-trainSize
trainIndex,testIndex = torch.utils.data.random_split(range(dataSize), [trainSize, testSize], generator=torch.Generator().manual_seed(2020))

xTrain = xss[trainIndex,:]
yTrain = yss[trainIndex,:]
xTest = xss[testIndex,:]
yTest = yss[testIndex,:]

# this is one way to define a network
class Net(torch.nn.Module):
    def __init__(self, n_feature, n_hidden, n_output):
        super(Net, self).__init__()
        self.hidden1 = torch.nn.Linear(n_feature, n_hidden[0])   # hidden layer
        self.hidden2 = torch.nn.Linear(n_hidden[0], n_hidden[1])   # hidden layer
        #self.hidden3 = torch.nn.Linear(n_hidden[1], n_hidden[2])   # hidden layer
        self.predict = torch.nn.Linear(n_hidden[1], n_output)       # output layer
        # Define proportion or neurons to dropout
        self.dropout = torch.nn.Dropout(0.8)
    def forward(self, x):
        x = torch.sigmoid(self.hidden1(x))      # activation function for hidden layer
        x = torch.sigmoid(self.hidden2(x))      # activation function for hidden layer
        #x = self.dropout(x)
        #x = self.dropout(x)
        #x = torch.sigmoid(self.hidden3(x))      # activation function for hidden layer
        #x = self.dropout(x)
        out = self.predict(x)               # linear output
        return out

# model architecture
hlayers = [30,10]
nfeatures = np.shape(xTrain)[1]
noutput = np.shape(yTrain)[1]
Arch = '_30_10'

if TRAIN:
	# torch can only train on Variable, so convert them to Variable
    xTrain = Variable(torch.from_numpy(xTrain))
    yTrain = Variable(torch.from_numpy(yTrain))
    
    # Check cuda availability
    cuda = torch.cuda.is_available()
    
    # Create neural network model
    if cuda:
        torch.cuda.manual_seed(2020)
        model = Net(n_feature = nfeatures, n_hidden = hlayers, n_output = noutput).cuda()
        device = 'cuda'
    else:
        torch.manual_seed(2020)
        model = Net(n_feature = nfeatures, n_hidden = hlayers, n_output = noutput)
        device = 'cpu'

    optimizer = torch.optim.LBFGS(model.parameters())
    loss_func = torch.nn.MSELoss()  # this is for regression mean squared loss
    
    model = model.double()
    xTrain = xTrain.to(device) 
    yTrain = yTrain.to(device)
        
    def closure():
        optimizer.zero_grad()
        output = model(xTrain)
        loss = loss_func(output, yTrain)
        loss.backward()
        return loss
    
    for epoch in range(200):                            # loop over the dataset multiple times
        prediction = model(xTrain)                      # input x and predict based on x

        #loss = loss_func(prediction, yTrain)     # must be (1. nn output, 2. target)

        #optimizer.zero_grad()          # clear gradients for next train
        #loss.backward()                # backpropagation, compute gradients
        loss = optimizer.step(closure)        # apply gradients

        # print statistics
        running_loss = loss.item()
        print('[%d] loss: %.7f' %(epoch + 1,running_loss))

    print('Finished Training')
	
    torch.save(model.state_dict(),'model%s.pth' %(Arch))

net = Net(n_feature=np.shape(x)[1], n_hidden=hlayers, n_output=np.shape(ati)[1])
net = net.double()
net.load_state_dict(torch.load('model%s.pth' %(Arch)))
pred = net(Variable(torch.from_numpy(x)))
pred = pred.detach().numpy()
pred = pred*max_ati

fig, axs = plt.subplots(3, 2)
axs[0, 0].plot(ati[:,0])
axs[0, 0].plot(pred[:,0])

axs[1, 0].plot(ati[:,1])
axs[1, 0].plot(pred[:,1])

axs[2, 0].plot(ati[:,2])
axs[2, 0].plot(pred[:,2])

axs[0, 1].plot(ati[:,3])
axs[0, 1].plot(pred[:,3])

axs[1, 1].plot(ati[:,4])
axs[1, 1].plot(pred[:,4])

axs[2, 1].plot(ati[:,5])
axs[2, 1].plot(pred[:,5])
plt.show()
