import torch


class LinearLayer(torch.nn.Module):
    '''A linear transformation.'''

    def __init__(self, input_size):
        super().__init__()
        self.linear = torch.nn.Linear(input_size, 1)
        self.input_size = input_size

    def forward(self, x):
        return self.linear(x)

class NN(torch.nn.Module):

    def __init__(self, input_size):
        super().__init__()
        self.linear = torch.nn.Linear(input_size, 1)
        self.input_size = input_size

    def forward(self, x):
        return torch.tanh(self.linear(x))
