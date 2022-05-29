# This code loads a trained network and tests it

import os
import torch
import torch.nn as nn
import numpy as np
import h5py
import time

# Load Data
#from google.colab import drive
#drive.mount('/content/drive')


#--------------- USER INPUT -------------------#
model_path = "./networks/NN1.pt"
#----------------------------------------------#



#------- Load model --------#
# initialize the classes
  # Standard Layer Definition
class StandardLayer(nn.Module):
  def __init__(self, n_input_features, n_output_features):
    super(StandardLayer, self).__init__()
    self.lin1 = nn.Linear(n_input_features, n_input_features)
    self.lin2 = nn.Linear(n_input_features, n_output_features)
    self.LRelu = nn.LeakyReLU(0.01)

  def forward(self, x):
    z1 = self.LRelu(self.lin1((x)))
    z2 = self.LRelu(self.lin2(z1)+self.lin2(x))
    return z2

  
Model = nn.Sequential(   # Layer 0
                      StandardLayer(55,40),  # Layer 1
                      StandardLayer(40,40),  # Layer 5
                      StandardLayer(40,40),  # Layer 5
                      StandardLayer(40,20),  # Layer 6
                      StandardLayer(20,20),  # Layer 10
                      StandardLayer(20,20),  # Layer 10
                      StandardLayer(20,10),  # Layer 12
                      StandardLayer(10,10),  # Layer 15
                      StandardLayer(10,10),  # Layer 15
                      StandardLayer(10,10),  # Layer 15
                      nn.Linear(10, 2))      # Layer 11 - Output Layer
Model.double()

#Model = nn.Sequential() 
#Model.double()
#Optimizer = torch.optim.Adam(Model.parameters())

# Read model from path
# Note this way of saving and loading will allow you 
# to continue training the model if required.
#Model = TheModelClass(*args, **kwargs)
#Optimizer = TheOptimizerClass(*args, **kwargs)

checkpoint = torch.load(model_path)
Model.load_state_dict(checkpoint['model_state_dict'])
#Optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
#epoch = checkpoint['epoch']

#Remember that you must call model.eval() to set dropout and batch normalization layers to evaluation mode before running inference. Failing to do this will yield inconsistent inference results. If you wish to resuming training, call model.train() to ensure these layers are in training mode.
Model.eval()
# - or -
#model.train()



#------------ evaluating errors -------#

# Error Evaluation of Model on Dev Data
print('evaluating model on dev data')
num_dev_batches = 900
#dev_dir = "/content/drive/MyDrive/training_data/test_selection_2/"
#dev_dir = "../../training_data/test_selection_2_onebatch/"
dev_dir = "../../training_data/test_selection_2/"
with torch.no_grad():
   norm_error = 0
   tau_error = 0
   q_error = 0
   for j in range(1, num_dev_batches+1):
      # Load test data from batch (i) into numpy here
      filename = dev_dir+"batch_"+ str(j).zfill(5) + ".hdf"
      f = h5py.File(filename, 'r')
      data = np.transpose(f['data'])
      Y_tau_np = data[:,0]
      Y_q_np = data[:,1]
      X_np = data[:,2:57]
      X = torch.from_numpy(X_np)
      Y_tau = torch.from_numpy(Y_tau_np)
      Y_q = torch.from_numpy(Y_q_np)
      Y_hat = Model(X)
      Y_tau_hat = Y_hat[:,0]
      Y_q_hat = Y_hat[:,1]
      q_error += 2*torch.mean(torch.div(torch.abs(Y_q-Y_q_hat),torch.add(torch.abs(Y_q),torch.abs(Y_q_hat))))
      tau_error += 2*torch.mean(torch.div(torch.abs(Y_tau-Y_tau_hat),torch.add(torch.abs(Y_tau),torch.abs(Y_tau_hat))))
   tau_error /= num_dev_batches
   q_error /= num_dev_batches
print("Dev q error: ", q_error.item())
print("Dev tau error: ", tau_error.item())





# Error Evaluation of Model on Test Data
print('evaluating model on test data')
#test_dir = "/content/drive/MyDrive/training_data/test_selection_2/"
#test_dir = "../../training_data/test_selection_2_onebatch/"
test_dir = "../../training_data/test_selection_2/"
with torch.no_grad():
   tau_error = 0
   q_error = 0
   norm_error = 0
   for j in range(num_dev_batches+1, 2*num_dev_batches+1):
      # Load test data from batch (i) into numpy here
      filename = test_dir+"batch_"+ str(j).zfill(5) + ".hdf"
      f = h5py.File(filename, 'r')
      data = np.transpose(f['data'])
      Y_tau_np = data[:,0]
      Y_q_np = data[:,1]
      X_np = data[:,2:57]
      X = torch.from_numpy(X_np)
      Y_tau = torch.from_numpy(Y_tau_np)
      Y_q = torch.from_numpy(Y_q_np)
      Y_hat = Model(X)
      Y_tau_hat = Y_hat[:,0]
      Y_q_hat = Y_hat[:,1]
      q_error += 2*torch.mean(torch.div(torch.abs(Y_q-Y_q_hat),torch.add(torch.abs(Y_q),torch.abs(Y_q_hat))))
      tau_error += 2*torch.mean(torch.div(torch.abs(Y_tau-Y_tau_hat),torch.add(torch.abs(Y_tau),torch.abs(Y_tau_hat))))
   tau_error /= num_dev_batches
   q_error /= num_dev_batches
print("Test q error: ", q_error.item())
print("Test tau error: ", tau_error.item())



