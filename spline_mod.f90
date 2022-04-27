module spline_mod
    use nrtype 
    implicit none 

    private:: tridag_ser,tridag,locate
    public:: spline,splint

    INTERFACE tridag_ser
        module procedure tridag_ser
    END INTERFACE

    INTERFACE tridag
        module procedure tridag_par
    END INTERFACE

    INTERFACE locate
        module procedure locate
    END INTERFACE

    INTERFACE spline
        module procedure spline
    END INTERFACE

    INTERFACE splint
        module procedure splint
    END INTERFACE


contains 

    SUBROUTINE tridag_ser(a,b,c,r,u)
    USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: a,b,c,r
    REAL(SP), DIMENSION(:), INTENT(OUT) :: u
    REAL(SP), DIMENSION(size(b)) :: gam
    INTEGER(I4B) :: n,j
    REAL(SP) :: bet
    n=assert_eq((/size(a)+1,size(b),size(c)+1,size(r),size(u)/),'tridag_ser')
    bet=b(1)
    if (bet == 0.0) call nrerror('tridag_ser: Error at code stage 1')
    u(1)=r(1)/bet
    do j=2,n
		gam(j)=c(j-1)/bet
		bet=b(j)-a(j-1)*gam(j)
		if (bet == 0.0) &
			call nrerror('tridag_ser: Error at code stage 2')
		u(j)=(r(j)-a(j-1)*u(j-1))/bet
	end do
	do j=n-1,1,-1
		u(j)=u(j)-gam(j+1)*u(j+1)
	end do
	END SUBROUTINE tridag_ser

	RECURSIVE SUBROUTINE tridag_par(a,b,c,r,u)
	USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
	!USE nr, ONLY : tridag_ser
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: a,b,c,r
	REAL(SP), DIMENSION(:), INTENT(OUT) :: u
	INTEGER(I4B), PARAMETER :: NPAR_TRIDAG=4
	INTEGER(I4B) :: n,n2,nm,nx
	REAL(SP), DIMENSION(size(b)/2) :: y,q,piva
	REAL(SP), DIMENSION(size(b)/2-1) :: x,z
	REAL(SP), DIMENSION(size(a)/2) :: pivc
	n=assert_eq((/size(a)+1,size(b),size(c)+1,size(r),size(u)/),'tridag_par')
	if (n < NPAR_TRIDAG) then
		call tridag_ser(a,b,c,r,u)
	else
		if (maxval(abs(b(1:n))) == 0.0) &
			call nrerror('tridag_par: possible singular matrix')
		n2=size(y)
		nm=size(pivc)
		nx=size(x)
		piva = a(1:n-1:2)/b(1:n-1:2)
		pivc = c(2:n-1:2)/b(3:n:2)
		y(1:nm) = b(2:n-1:2)-piva(1:nm)*c(1:n-2:2)-pivc*a(2:n-1:2)
		q(1:nm) = r(2:n-1:2)-piva(1:nm)*r(1:n-2:2)-pivc*r(3:n:2)
		if (nm < n2) then
			y(n2) = b(n)-piva(n2)*c(n-1)
			q(n2) = r(n)-piva(n2)*r(n-1)
		end if
		x = -piva(2:n2)*a(2:n-2:2)
		z = -pivc(1:nx)*c(3:n-1:2)
		call tridag_par(x,y,z,q,u(2:n:2))
		u(1) = (r(1)-c(1)*u(2))/b(1)
		u(3:n-1:2) = (r(3:n-1:2)-a(2:n-2:2)*u(2:n-2:2) &
			-c(3:n-1:2)*u(4:n:2))/b(3:n-1:2)
		if (nm == n2) u(n)=(r(n)-a(n-1)*u(n-1))/b(n)
	end if
	END SUBROUTINE tridag_par


    FUNCTION locate(xx,x)
	USE nrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: xx
	REAL(SP), INTENT(IN) :: x
	INTEGER(I4B) :: locate
	INTEGER(I4B) :: n,jl,jm,ju
	LOGICAL :: ascnd
	n=size(xx)
	ascnd = (xx(n) >= xx(1))
	jl=0
	ju=n+1
	do
		if (ju-jl <= 1) exit
		jm=(ju+jl)/2
		if (ascnd .eqv. (x >= xx(jm))) then
			jl=jm
		else
			ju=jm
		end if
	end do
	if (x == xx(1)) then
		locate=1
	else if (x == xx(n)) then
		locate=n-1
	else
		locate=jl
	end if
	END FUNCTION locate


    SUBROUTINE spline(x,y,yp1,ypn,y2)
    USE nrtype; USE nrutil, ONLY : assert_eq
    !USE nr, ONLY : tridag
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
    REAL(SP), INTENT(IN) :: yp1,ypn
    REAL(SP), DIMENSION(:), INTENT(OUT) :: y2
    INTEGER(I4B) :: n
    REAL(SP), DIMENSION(size(x)) :: a,b,c,r
    n=assert_eq(size(x),size(y),size(y2),'spline')
    c(1:n-1)=x(2:n)-x(1:n-1)
    r(1:n-1)=6.0_sp*((y(2:n)-y(1:n-1))/c(1:n-1))
    r(2:n-1)=r(2:n-1)-r(1:n-2)
    a(2:n-1)=c(1:n-2)
    b(2:n-1)=2.0_sp*(c(2:n-1)+a(2:n-1))
    b(1)=1.0
    b(n)=1.0
    
    if (yp1 > 0.99e30_sp) then
        r(1)=0.0
        c(1)=0.0
    else
        r(1)=(3.0_sp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
        c(1)=0.5
    end if
    
    if (ypn > 0.99e30_sp) then
        r(n)=0.0
        a(n)=0.0
    else
        r(n)=(-3.0_sp/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)
        a(n)=0.5
    end if

    call tridag(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))   
    END SUBROUTINE spline
    
    FUNCTION splint(xa,ya,y2a,x)
    USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
    !USE nr, ONLY: locate
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: xa,ya,y2a
    REAL(SP), INTENT(IN) :: x
    REAL(SP) :: splint
    INTEGER(I4B) :: khi,klo,n
    REAL(SP) :: a,b,h
    n=assert_eq(size(xa),size(ya),size(y2a),'splint')
    klo=max(min(locate(xa,x),n-1),1)
    khi=klo+1
    h=xa(khi)-xa(klo)
    if (h == 0.0) call nrerror('bad xa input in splint')
    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    splint=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0_sp
    END FUNCTION splint

END MODULE spline_mod