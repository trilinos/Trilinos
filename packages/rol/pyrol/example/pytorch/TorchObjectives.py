from pyrol import Objective

import torch


class TorchObjective(Objective):
    # https://pytorch.org/docs/stable/func.html

    # @staticmethod
    # def _copy(source, target):
    #     target.zero()
    #     target.plus(source)

    def __init__(self):
        super().__init__()
        self.torch_gradient = torch.func.grad(self.torch_value)

    def value(self, x, tol):
        return self.torch_value(x.torch_object).item()

    def gradient(self, g, x, tol):
        ans = self.torch_gradient(x.torch_object)
        g.torch_object = ans

    def _forward_over_reverse(self, input, x, v):
        # https://github.com/google/jax/blob/main/docs/notebooks/autodiff_cookbook.ipynb
        return torch.func.jvp(input, (x,), (v,))

    def hessVec(self, hv, v, x, tol):
        input = torch.func.grad(self.torch_value)
        _, ans = self._forward_over_reverse(input, x.torch_object, v.torch_object)
        hv.torch_object = ans

    def torch_value(self, x):
        # Returns a scalar torch Tensor
        raise NotImplementedError


class SquaredNorm(TorchObjective):

    def torch_value(self, x):
        return 0.5*x.squeeze()**2


class LeastSquaresObjective(TorchObjective):

    def __init__(self, data, model):
        super().__init__()
        self.x, self.y = data
        self.model = model
        self.loss = torch.nn.MSELoss(reduction='sum')

    def torch_value(self, x):
        return 0.5*self.loss(torch.func.functional_call(self.model, x, self.x), self.y)
