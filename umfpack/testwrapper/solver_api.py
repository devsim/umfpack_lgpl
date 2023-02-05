from enum import Enum
import abc

class Format(Enum):
    CCM = 1
    CRM = 2

class Preconditioner:
    def __init__(self, numrows, matrix_format):
        self.numrows = numrows
        self.matrix_format = matrix_format


class UMFPACK(Preconditioner):
    def __init__(self, numrows):
        super(UMFPACK, self).__init__(numrows, Format.CCM)

foo = UMFPACK(5)
print(foo.matrix_format)
