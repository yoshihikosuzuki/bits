from abc import ABC, abstractmethod
from typing import List


class ExplicitRepr(ABC):
    """An abstract class where the order of variables displayed in the repr()
    string must be explicityly specified.

    The `__repr__` method must be defined in every child class to specify the
    order of variables to be displayed in the repr() string.
    Pass a list of variable names to the helper method, `_order_repr()`, to
    generate the repr() string.
    """

    def _order_repr(self, var_names: List[str]) -> str:
        var_reprs = ', '.join(map(lambda x: f"{x}={repr(getattr(self, x))}",
                                  var_names))
        return f"{self.__class__.__name__}({var_reprs})"

    @abstractmethod
    def __repr__(self) -> str:
        return self._order_repr(["seq"])   # NOTE: This is an example
