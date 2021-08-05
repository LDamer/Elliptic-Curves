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

    def multiply_with_scalar(self, p, e):
        l = len(bin(e)[2:])
        res = (self.n, self.n)
        for i in range(l):
            res = self.add(res, res)
            if e & (1 << (l - i - 1)) == 1 << (l - i - 1):
                res = self.add(p, res)
        return res

    def add(self, p, q):
        """
        :param p: point on the curve. Negative values allowed.
        :param q: point on the curve. Negative values allowed.
        """
        # Check for correct inputs. Both inputs must be of type tupel.
        if type(p) != type(q) or type(p) != tuple:
            raise Exception("x and y must be tuples!")

        # p must be part of the curve or has to be the neutral element.
        # In case it is the neutral element we return q
        if p != (self.n, self.n):
            if (p[0] ** 3 + self.a * p[0] + self.b) % self.p != (p[1] ** 2) % self.p:
                raise Exception("p is not part of the curve!")
        else:
            return q

        # q must be part of the curve or has to be the neutral element.
        # In case it is the neutral element we return p
        if q != (self.n, self.n):
            if (q[0] ** 3 + self.a * q[0] + self.b) % self.p != (q[1] ** 2) % self.p:
                raise Exception("q is not part of the curve!")
        else:
            return p

        # Check whether the points are the inverse of each other.
        # Every comparison just after modulo reduction for better usability (negative values allowed).
        # Therefore, you can directly type the inverse of P as (p.x, -p.y) without reducing by yourself modulo p.
        if p[0] % self.p == q[0] % self.p and (-p[1]) % self.p == q[1] % self.p:
            return self.n, self.n
        else:
            # Now it is guarenteed that neither of the points is the neutral Element and they both lie on the curve.
            if p[0] % self.p == q[0] % self.p and p[1] % self.p == q[1] % self.p:
                # In this case we need to perform point multiplication.
                # A lot of modulo reduction because negative values are allowed.
                s = ((3 * p[0] ** 2 + self.a) * pow(2 * p[1] % self.p, -1, self.p)) % self.p
                x_3 = (s ** 2 - p[0] - q[0]) % self.p
                y_3 = (s * (p[0] - x_3) - p[1]) % self.p
                return x_3, y_3
            else:
                # In this case we need to perform point addition.
                # A lot of modulo reduction because negative values are allowed.
                s = (((q[1] - p[1]) % self.p) * pow((q[0] % self.p) - (p[0] % self.p), -1, self.p)) % self.p
                x_3 = (s ** 2 - p[0] - q[0]) % self.p
                y_3 = (s * (p[0] - x_3) - p[1]) % self.p
                return x_3, y_3

