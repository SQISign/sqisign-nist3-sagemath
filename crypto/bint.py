from sage.rings.integer import Integer as _SageInteger

import operator

_CASTABLE_TYPES = (int, _SageInteger)

class BoundedInt:
    r"""
    An integer with a declared signed bitlen bound [Section 4.1.2].

    Bounds are propagated through arithmetic as in the reference
    implementation, so that sampling random numbers is done consistently with
    the KATs

    INPUT:

    - ``n`` -- an integer
    - ``b`` -- the bound (default: the bit length of `|n|` plus one)
    """

    def __init__(self, n, b=None):
        if isinstance(n, BoundedInt):
            return n

        if not b:
            b = abs(int(n)).bit_length() + 1
        self.bb = b
        self.nn = int(n)

    @classmethod
    def _cast(cls, other):
        if isinstance(other, BoundedInt):
            return other
        if isinstance(other, _CASTABLE_TYPES):
            return cls(other)
        raise NotImplementedError(f'{type(other)} not supported')

    def __int__(self):
        return self.nn

    def __repr__(self):
        return f'BoundedInt({self.nn}, {self.bb})'

    def __str__(self):
        return self.__repr__()

    def __add__(self, other):
        other = self._cast(other)
        n1 = self.nn + other.nn
        b1 = max(self.bb, other.bb) + 1
        return BoundedInt(n1, b1)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        other = self._cast(other)
        n1 = self.nn - other.nn
        b1 = max(self.bb, other.bb) + 1
        return BoundedInt(n1, b1)

    def __rsub__(self, other):
        other = self._cast(other)
        n1 = other.nn - self.nn
        b1 = max(self.bb, other.bb) + 1
        return BoundedInt(n1, b1)

    def __neg__(self):
        return BoundedInt(-self.nn, self.bb)

    def __mul__(self, other):
        other = self._cast(other)
        n1 = self.nn * other.nn
        b1 = self.bb + other.bb - 1
        return BoundedInt(n1, b1)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __hash__(self):
        return hash(self.nn)

    def _compare(self, other, op):
        other = self._cast(other)
        return op(self.nn, other.nn)

    def __eq__(self, other):
        return self._compare(other, operator.eq)

    def __ne__(self, other):
        return self._compare(other, operator.ne)

    def __lt__(self, other):
        return self._compare(other, operator.lt)

    def __le__(self, other):
        return self._compare(other, operator.le)

    def __gt__(self, other):
        return self._compare(other, operator.gt)

    def __ge__(self, other):
        return self._compare(other, operator.ge)

    def _integer_(self, ZZ):
        return ZZ(self.nn)
