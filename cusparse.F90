module cusparse

interface Dgtsv
         module procedure cusparseDtgsv
         module procedure cusparseStgsv
end interface Dgtsv


interface cusparseDtgsv
subroutine cusparseDtgsv()
end subroutine
end interface

interface cusparseStgsv
subroutine cusparseStgsv()
end subroutine
end interface


end module cusparse
