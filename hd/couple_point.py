from ec.ec_point import ECPoint

class CoupleECPoint:
    r"""
    Implements the type ECPointPair [Definition 4.5.6].

    A point `(P_1, P_2)` on a product `E_1 \times E_2` of Montgomery curves.

    INPUT:

    - ``P1`` -- a point on `E_1`
    - ``P2`` -- a point on `E_2`
    """

    def __init__(self, P1: ECPoint, P2: ECPoint):
        self.P1 = P1
        self.P2 = P2

    def __repr__(self):
        return "[{}, {}]".format(self.P1, self.P2)

    def points(self):
        return self.P1, self.P2

    def dbl(self):
        P1, P2 = self.points()
        P1 = P1.ec_dbl()
        P2 = P2.ec_dbl()
        return CoupleECPoint(P1, P2)
