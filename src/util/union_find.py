from dataclasses import dataclass, field
from typing import List


@dataclass(eq=False)
class UnionFind:
    N: int
    parents: List[int] = field(init=None)

    def __post_init__(self):
        self.parents = list(range(self.N))

    def get_root(self, x):
        if self.parents[x] == x:
            return x
        self.parents[x] = self.get_root(self.parents[x])
        return self.parents[x]

    def in_same_set(self, x, y):
        return self.get_root(x) == self.get_root(y)

    def unite(self, x, y):
        x = self.get_root(x)
        y = self.get_root(y)
        if x == y:
            return
        self.parents[x] = y
