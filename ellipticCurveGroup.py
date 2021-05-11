class ECurve:
    def __init__(self, a, b, p, n):
        """
        a, b and p are the arguments for the elliptic curve of the form: y^2 = x^3 + a * x^2 + b mod p.  
        n represents the value of the neutral element. For userfriendly reasons n has to be a number less than 0.
        The neutral element itself is then represented as the tupel (n, n).
        Therefore, in usage the neutral must be given as tupel and not just as the integer value n! 
        e.g: n = -1 --> ECurve.add( (x_1, y_1), (-1,-1) ) returns (x_1, y_1) as result.
        """
        if 0 <= n:
            raise Exception("neutral element must be negative!")
        elif (4 * a ** 3 + 27 * b ** 2) % p == 0:
            raise Exception("this elliptic curve has singularities! Discriminant condition failed!")
        else:
            self.n = n
            self.p = p
            self.b = b
            self.a = a

    def ord_group(self):
        y = []
        counter = 0
        for i in range(self.p):
            y.append((i ** 2) % self.p)
        for x in range(self.p):
            y_s = (x ** 3 + self.a * x + self.b) % self.p
            for j in range(len(y)):
                tmp = 0
                if y_s == y[j]:
                    counter += 1
                    tmp += 1
                    if tmp == 2:
                        break
        return counter + 1

    def ord_element(self, a):
        result = a
        counter = 1
        while not result == (self.n, self.n):
            result = self.add(result, a)
            counter += 1
        return counter

    def add(self, P, Q):
        # Adjust the values of P and Q to the underliying finite field
        # in case one of their coordinates is negative
        P[0] = P[0] % self.p
        P[0] = P[1] % self.p
        Q[0] = Q[0] % self.p
        Q[0] = Q[1] % self.p

        # Check for correct inputs. Both inputs must be of type tupel.
        if type(P) != type(Q) or type(P) != tuple:
            raise Exception("x and y must be tuples!")

        # P must be part of the curve or has to be the neutral element.
        if P != (self.n, self.n):
            if (P[0] ** 3 + self.a * P[0] + self.b) % self.p != (P[1] ** 2) % self.p:
                raise Exception("P is not part of the curve!")

        # Q must be part of the curve or has to be the neutral element.
        if Q != (self.n, self.n):
            if (Q[0] ** 3 + self.a * Q[0] + self.b) % self.p != (Q[1] ** 2) % self.p:
                raise Exception("Q is not part of the curve!")

        # Check whether P is the neutral element.
        if P == (self.n, self.n):
            return Q

        # Check whether Q is the neutral element.
        elif Q == (self.n, self.n):
            return P

        # Check whether the points are the inverse of each other.
        elif P[0] == Q[0] and -P[1] % self.p == Q[1]:
            return self.n, self.n
        else:
            # Now is guarenteed that neither of the points is the neutral Element and they both lie on the curve.
            if P == Q:
                # In this case we need to perform point multiplication.
                s = ((3 * P[0] ** 2 + self.a) * pow(2 * P[1], -1, self.p)) % self.p
                x_3 = (s ** 2 - P[0] - Q[0]) % self.p
                y_3 = (s * (P[0] - x_3) - P[1]) % self.p
                return x_3, y_3
            else:
                # In this case we nee to perform point addition.               
                s = (((Q[1] - P[1]) % self.p) * pow(Q[0] - P[0], -1, self.p)) % self.p
                x_3 = (s ** 2 - P[0] - Q[0]) % self.p
                y_3 = (s * (P[0] - x_3) - P[1]) % self.p
                return x_3, y_3
