from pyrol.pyrol import ROL
from pyrol.getTypeName import *

import torch

import copy


class PythonVector(getTypeName('Vector')):

    def __iadd__(self, other):
        self.plus(other)

    # TO-DO: Add other operations.


class TensorVector(PythonVector):

    @torch.no_grad()
    def __init__(self, tensor):
        super().__init__()
        assert isinstance(tensor, torch.Tensor)
        self.torch_object = tensor

    @property
    def tensor(self):
        return self.torch_object

    @torch.no_grad()
    def axpy(self, alpha, other):
        self.tensor.add_(other.tensor, alpha=alpha)

    @torch.no_grad()
    def scale(self, alpha):
        self.tensor.mul_(alpha)

    @torch.no_grad()
    def zero(self):
        self.tensor.zero_()

    @torch.no_grad()
    def dot(self, other):
        ans = torch.sum(torch.mul(self.tensor, other.tensor))
        return ans.item()

    @torch.no_grad()
    def clone(self):
        tensor = copy.deepcopy(self.tensor)
        ans = TensorVector(tensor)
        ans.zero()
        return ans

    @torch.no_grad()
    def dimension(self):
        return self.tensor.numel()

    @torch.no_grad()
    def setScalar(self, alpha):
        self.fill_(alpha)

    @torch.no_grad()
    def __getitem__(self, index):
        flat = self.tensor.view(-1)
        return flat[index].item()

    @torch.no_grad()
    def __setitem__(self, index, value):
        flat = self.tensor.view(-1)
        flat[index] = value

    # Derived methods #########################################################

    @torch.no_grad()
    def reduce(self, op):
        reduction_type = op.reductionType()
        match reduction_type:
            case ROL.Elementwise.REDUCE_MIN:
                ans = self.tensor.min()
                ans = ans.item()
            case ROL.Elementwise.REDUCE_MAX:
                ans = self.tensor.max()
                ans = ans.item()
            case ROL.Elementwise.REDUCE_SUM:
                ans = torch.sum(self.tensor)
                ans = ans.item()
            case ROL.Elementwise.REDUCE_AND:
                ans = self.tensor.all()
                ans = ans.item()
            case ROL.Elementwise.REDUCE_BOR:
                ans = 0
                for i in range(self.dimension()):
                    ans = ans | int(self[i].item())
            case _:
                raise NotImplementedError(reduction_type)
        return ans

    @torch.no_grad()
    def plus(self, other):
        self.axpy(1, other)

    @torch.no_grad()
    def norm(self):
        return self.dot(self)**0.5

    @torch.no_grad()
    def basis(self, i):
        b = self.clone()
        b.zero()
        b[i] = 1
        return b

    ####

    @torch.no_grad()
    def applyUnary(self, op):
        for i in range(self.dimension()):
            self[i] = op.apply(self[i])

    @torch.no_grad()
    def applyBinary(self, other, op):
        for i in range(self.dimension()):
            self[i] = op.apply(self[i], other[i])


class TensorDictVector(PythonVector):

    @torch.no_grad()
    def __init__(self, tensor_dict):
        super().__init__()
        assert isinstance(tensor_dict, dict)
        self.torch_object = tensor_dict

    @property
    def tensor_dict(self):
        return self.torch_object

    @torch.no_grad()
    def axpy(self, alpha, other):
        for k, v in self.tensor_dict.items():
            v.add_(other.tensor_dict[k], alpha=alpha)

    @torch.no_grad()
    def scale(self, alpha):
        for _, v in self.tensor_dict.items():
            v.mul_(alpha)

    @torch.no_grad()
    def zero(self):
        for _, v in self.tensor_dict.items():
            v.zero_()

    @torch.no_grad()
    def dot(self, other):
        ans = 0
        for k, v in self.tensor_dict.items():
            ans += torch.sum(torch.mul(v, other.tensor_dict[k]))
        return ans.item()

    @torch.no_grad()
    def clone(self):
        tensor_dict = copy.deepcopy(self.tensor_dict)
        ans = TensorDictVector(tensor_dict)
        ans.zero()
        return ans

    @torch.no_grad()
    def dimension(self):
        # TO-DO: Cache value
        ans = 0
        for _, v in self.tensor_dict.items():
            ans += v.numel()
        return ans

    @torch.no_grad()
    def setScalar(self, alpha):
        for _, v in self.tensor_dict.items():
            v.fill_(alpha)

    @torch.no_grad()
    def __getitem__(self, index):
        total = 0
        for _, v in self.tensor_dict.items():
            numel = v.numel()
            if index < numel + total:
                flat = v.view(-1)
                return flat[index - total].item()
            total += numel

    @torch.no_grad()
    def __setitem__(self, index, value):
        total = 0
        for _, v in self.tensor_dict.items():
            numel = v.numel()
            if index < numel + total:
                flat = v.view(-1)
                flat[index - total] = value
                return
            total += numel

    # Derived methods #########################################################

    @torch.no_grad()
    def reduce(self, op):
        reduction_type = op.reductionType()
        match reduction_type:
            case ROL.Elementwise.REDUCE_MIN:
                ans = float('+inf')
                for _, v in self.torch_dict.items():
                    ans = min(ans, v.min().item())
            case ROL.Elementwise.REDUCE_MAX:
                ans = float('-inf')
                for _, v in self.torch_dict.items():
                    ans = max(ans, v.max().item())
            case ROL.Elementwise.REDUCE_SUM:
                ans = 0
                for _, v in self.torch_dict.items():
                    ans += torch.sum(v)
                ans = ans.item()
            case ROL.Elementwise.REDUCE_AND:
                ans = True
                for _, v in self.torch_dict.items():
                    ans = ans and v.all().item()
                    if ans == False:
                        break
            case ROL.Elementwise.REDUCE_BOR:
                ans = 0
                for i in range(self.dimension()):
                    ans = ans | int(self[i].item())
            case _:
                raise NotImplementedError(reduction_type)
        return ans

    @torch.no_grad()
    def plus(self, other):
        self.axpy(1, other)

    @torch.no_grad()
    def norm(self):
        return self.dot(self)**0.5

    @torch.no_grad()
    def basis(self, i):
        b = self.clone()
        b.zero()
        b[i] = 1
        return b

    ####

    @torch.no_grad()
    def applyUnary(self, op):
        for i in range(self.dimension()):
            self[i] = op.apply(self[i])

    @torch.no_grad()
    def applyBinary(self, other, op):
        for i in range(self.dimension()):
            self[i] = op.apply(self[i], other[i])
