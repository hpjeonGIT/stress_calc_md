!
!############################################################################
SUBROUTINE CELL_INIT(Nset, cell)
USE DATAFMT
IMPLICIT NONE
!
TYPE(AMNT):: Nset
TYPE(SORT):: cell(Nset%Ncell_all)
!
INTEGER(GI):: i, j, k, Ncell(3), ID_CELL, id
Ncell(:) = Nset%Ncell(:)
!
DO i=1, Ncell(1)
   DO j=1, Ncell(2)
      DO k=1, Ncell(3)
         id = ID_CELL(Ncell, i, j, k)
         cell(id)%pair( 1) = ID_CELL(Ncell, i-1, j-1, k-1)
         cell(id)%pair( 2) = ID_CELL(Ncell, i  , j-1, k-1)
         cell(id)%pair( 3) = ID_CELL(Ncell, i+1, j-1, k-1)
         cell(id)%pair( 4) = ID_CELL(Ncell, i-1, j  , k-1)
         cell(id)%pair( 5) = ID_CELL(Ncell, i  , j  , k-1)
         cell(id)%pair( 6) = ID_CELL(Ncell, i+1, j  , k-1)
         cell(id)%pair( 7) = ID_CELL(Ncell, i-1, j+1, k-1)
         cell(id)%pair( 8) = ID_CELL(Ncell, i  , j+1, k-1)
         cell(id)%pair( 9) = ID_CELL(Ncell, i+1, j+1, k-1)
         cell(id)%pair(10) = ID_CELL(Ncell, i-1, j-1, k  )
         cell(id)%pair(11) = ID_CELL(Ncell, i  , j-1, k  )
         cell(id)%pair(12) = ID_CELL(Ncell, i+1, j-1, k  )
         cell(id)%pair(13) = ID_CELL(Ncell, i-1, j  , k  )
      END DO
   END DO
END DO
!
RETURN
END SUBROUTINE CELL_INIT
!
!*****************************************************************************
FUNCTION ID_CELL(Ncell, i, j, k)
USE DATAFMT, ONLY:GI
IMPLICIT NONE
!
INTEGER(GI):: Ncell(3), i, j, k, ID_CELL
INTEGER(GI):: l, m, n
l = MOD(i-1+Ncell(1), Ncell(1))
m = MOD(j-1+Ncell(2), Ncell(2))
n = MOD(k-1+Ncell(3), Ncell(3))
ID_CELL = l + m*Ncell(1) + n*Ncell(1)*Ncell(2) + 1
RETURN
END FUNCTION ID_CELL
!
!############################################################################
SUBROUTINE CELL_SORT(Nset, cell, q)
USE DATAFMT
TYPE(AMNT):: Nset
TYPE(SORT):: cell(Nset%Ncell_all)
TYPE(PTCL):: q(Nset%Npt_all)
!
INTEGER(GI):: i, Npt(Nset%Ncell_all), Npt_all, Ncell(3), id(3), id_cell
REAL   (DP):: hbox(3), Lcell(3), xx(3)
!
Npt_all = Nset%Npt_all
hbox(:) = 0.5D0*Nset%box(:)
Ncell(:) = Nset%Ncell(:)
Lcell(:) = Nset%Lcell(:)
!
Npt(:) = 0
DO i=1, Npt_all
   xx(:) = q(i)%xx(:) + hbox(:)
   id(:) = INT(xx(:)/Lcell(:))
   id(:) = MOD(id(:) + Ncell(:), Ncell(:))
   id_cell = id(1) + id(2)*Ncell(1) + id(3)*Ncell(1)*Ncell(2) + 1
   Npt(id_cell) = Npt(id_cell) + 1
   cell(id_cell)%link(Npt(id_cell)) = i
   IF (Npt(id_cell) > Npt_max) STOP "=== Cell particle limit crashed !!! ==="
END DO
cell(:)%Npt = Npt(:)
!
RETURN
END SUBROUTINE CELL_SORT
!
!############################################################################
SUBROUTINE VOLUME_CONTROL(Nset, sys, q, tag)
USE DATAFMT
IMPLICIT NONE
!
TYPE(AMNT):: Nset
TYPE(STAT):: sys
TYPE(PTCL):: q(Nset%Npt_all)
LOGICAL   :: tag
!
INTEGER(GI):: Ncell_new(3), i
REAL   (DP):: Lcell_new(3), box(3), box_new(3), P_int(3), lambda, &
     & lambda3d(3), P_short, xx(3), dx(3)
!
box(1:3) = sys%box(1:3)
P_int(:) = 0.0D0
!
! Pressure calculations
DO i=1, Nset%Npt_all
   P_int(:) = P_int(:) + q(i)%xx(:)*q(i)%ff(:)
END DO
!
! Page. 47, Allen & Tildesley, \sum r_i f_i = - \sum r_ij f_ij
P_int(:) = (sys%mvx2(:) + P_int(:))/(3.D0*Nset%V)
P_short = (P_int(1) + P_int(2) + P_int(3))/3.D0
sys%press = P_short
!
! Volume update
IF (Nset%iso_volume) THEN
   lambda = (1.D0 - Nset%beta_P*(P_short - sys%P_given))**INVTHR
   box_new(:) = box(:)/lambda
   lambda3d(:) = lambda
ELSE
   lambda3d(:) = (1.D0 - Nset%beta_P*(P_int(:) - sys%P_given))**INVTHR
   box_new(:) = box(:)/lambda3d(:)
END IF
Ncell_new(:) = INT(box_new(:)/sys%Rcut)
Lcell_New(:) = box_new(:)/DBLE(Ncell_new(:))
DO i=1, Nset%Npt_all
   xx(:) = q(i)%xx(:)
   q(i)%xx(:) = q(i)%xx(:)/lambda3d(:)
   dx(:) = q(i)%xx(:) - xx(:)
   q(i)%xr(:) = q(i)%xr(:) + dx(:)
END DO
!
tag = .FALSE.
IF (Ncell_new(1) /= Nset%Ncell(1) .OR. Ncell_new(2) /= Nset%Ncell(2) &
     & .OR. Ncell_new(3) /= Nset%Ncell(3)) THEN
   IF (Ncell_new(1) < 3 .AND. Ncell_new(2) < 3 .AND. Ncell_new(3) < 3) THEN
      PRINT *, "=== NPT ERROR - the volume is too small now ===", box_new(:)
      STOP
   END IF
   tag = .TRUE.
   Nset%Ncell(:) = Ncell_new(:)
   Nset%Ncell_all = Ncell_new(1)*Ncell_new(2)*Ncell_new(3)
END IF
Nset%box(:) = box_new(:)
Nset%Lcell(:) = Lcell_new(:)
sys%box(:) = box_new(:)
Nset%V = box_new(1)*box_new(2)*box_new(3)
!
END SUBROUTINE VOLUME_CONTROL
!
!############################################################################
SUBROUTINE HEINZ_STRESS(Nset, sys, param, q)
USE DATAFMT
IMPLICIT NONE
!
TYPE(AMNT):: Nset
TYPE(STAT):: sys
TYPE(ATMC):: param(Nset%Natom)
TYPE(PTCL):: q(Nset%Npt_all)
!
INTEGER(GI):: i, j, id, nd(3), md(3), N_grid(3)
REAL   (DP):: dx(3), dA(3), dA2(3), hbox(3), dv(3), xm, dL(3), r, r2, &
     & ss(sys%N_grid(1),sys%N_grid(2),sys%N_grid(3),9), d, d2

N_grid(:) = sys%N_grid(:)
dx(:) = sys%dx(:)
dA(:) = sys%dA(:)
hbox(:) = sys%box(:)*0.5D0
ss(:,:,:,:) = 0.0D0
dv(:) = 2.D0*dx(:)*dA(:)
dA2(:) = dA(:)*2.D0
DO i=1, Nset%Npt_all
   nd(:) = INT((q(i)%xx(:) + hbox(:))/dx(:)) + 1
   DO j=1,3
      IF (md(j) > N_grid(j)) md(j) = 1
   END DO
   id = q(i)%id
   xm = param(id)%am   
   !
   ! kinetic
   ss(nd(1),nd(2),nd(3),1) = ss(nd(1),nd(2),nd(3),1) - &
        &  xm*q(i)%xv(1)*q(i)%xv(1)/dv(1)
   ss(nd(1),nd(2),nd(3),2) = ss(nd(1),nd(2),nd(3),2) - &
        &xm*q(i)%xv(2)*q(i)%xv(2)/dv(2)
   ss(nd(1),nd(2),nd(3),3) = ss(nd(1),nd(2),nd(3),3) - &
        & xm*q(i)%xv(3)*q(i)%xv(3)/dv(3)
   ss(nd(1),nd(2),nd(3),4) = ss(nd(1),nd(2),nd(3),4) - &
        & xm*q(i)%xv(1)*q(i)%xv(2)/dv(2)
   ss(nd(1),nd(2),nd(3),5) = ss(nd(1),nd(2),nd(3),5) - &
        & xm*q(i)%xv(2)*q(i)%xv(3)/dv(3)
   ss(nd(1),nd(2),nd(3),6) = ss(nd(1),nd(2),nd(3),6) - &
        & xm*q(i)%xv(3)*q(i)%xv(1)/dv(1)
   !
   ss(md(1),md(2),md(3),1) = ss(md(1),md(2),md(3),1) - &
        &  xm*q(i)%xv(1)*q(i)%xv(1)/dv(1)
   ss(md(1),md(2),md(3),2) = ss(md(1),md(2),md(3),2) - &
        &xm*q(i)%xv(2)*q(i)%xv(2)/dv(2)
   ss(md(1),md(2),md(3),3) = ss(md(1),md(2),md(3),3) - &
        & xm*q(i)%xv(3)*q(i)%xv(3)/dv(3)
   ss(md(1),md(2),md(3),4) = ss(md(1),md(2),md(3),4) - &
        & xm*q(i)%xv(1)*q(i)%xv(2)/dv(2)
   ss(md(1),md(2),md(3),5) = ss(md(1),md(2),md(3),5) - &
        & xm*q(i)%xv(2)*q(i)%xv(3)/dv(3)
   ss(md(1),md(2),md(3),6) = ss(md(1),md(2),md(3),6) - &
        & xm*q(i)%xv(3)*q(i)%xv(1)/dv(1)
   !
   ! Force
   ss(nd(1),nd(2),nd(3),1) = ss(nd(1),nd(2),nd(3),1) + &
        &  q(i)%ff(1)/dA2(1)
   ss(nd(1),nd(2),nd(3),2) = ss(nd(1),nd(2),nd(3),2) + &
        &  q(i)%ff(2)/dA2(2)
   ss(nd(1),nd(2),nd(3),3) = ss(nd(1),nd(2),nd(3),3) + &
        &  q(i)%ff(3)/dA2(3)
   ss(nd(1),nd(2),nd(3),4) = ss(nd(1),nd(2),nd(3),4) + &
        &  q(i)%ff(1)/dA2(2)
   ss(nd(1),nd(2),nd(3),5) = ss(nd(1),nd(2),nd(3),5) + &
        &  q(i)%ff(2)/dA2(3)
   ss(nd(1),nd(2),nd(3),6) = ss(nd(1),nd(2),nd(3),6) + &
        &  q(i)%ff(3)/dA2(1)
   !
   ss(md(1),md(2),md(3),1) = ss(md(1),md(2),md(3),1) - &
        &  q(i)%ff(1)/dA2(1)
   ss(md(1),md(2),md(3),2) = ss(md(1),md(2),md(3),2) - &
        &  q(i)%ff(2)/dA2(2)
   ss(md(1),md(2),md(3),3) = ss(md(1),md(2),md(3),3) - &
        &  q(i)%ff(3)/dA2(3)
   ss(md(1),md(2),md(3),4) = ss(md(1),md(2),md(3),4) - &
        &  q(i)%ff(1)/dA2(2)
   ss(md(1),md(2),md(3),5) = ss(md(1),md(2),md(3),5) - &
        &  q(i)%ff(2)/dA2(3)
   ss(md(1),md(2),md(3),6) = ss(md(1),md(2),md(3),6) - &
        &  q(i)%ff(3)/dA2(1)
   !
   ! Radial and hoop stress
   dL(:) = hbox(:)
   DO j=1, 3
      IF (q(i)%xx(j) < 0.) dL(j) = -dL(j)
   END DO
   dL(:) = q(i)%xx(:) - dL(:)
   r2 = dL(1)*dL(1) + dL(2)*dL(2) + dL(3)*dL(3)
   d2 = dL(1)*dL(1) + dL(2)*dL(2)
   IF (r2 > 0.1 .AND. d2 > 0.1 ) THEN
      ! http://en.wikipedia.org/wiki/List_of_common_coordinate_transformations
      r = DSQRT(r2)
      d = DSQRT(d2)
      dL(:) = dL(:) / r
      ss(md(1),md(2),md(3),7) =  ss(md(1),md(2),md(3),7) + &
           & dL(1)*ss(md(1),md(2),md(3),1) + dL(2)*ss(md(1),md(2),md(3),2) + &
           & dL(3)*ss(md(1),md(2),md(3),3)
      ss(md(1),md(2),md(3),8) =  ss(md(1),md(2),md(3),8) + &
           & dL(1)*dL(3)*ss(md(1),md(2),md(3),1)/d + &
           & dL(2)*dL(3)*ss(md(1),md(2),md(3),2)/d - &
           & d*ss(md(1),md(2),md(3),3)/r2
      ss(md(1),md(2),md(3),9) =  ss(md(1),md(2),md(3),9) - &
           & dL(2)*r*ss(md(1),md(2),md(3),1)/d2 + &
           & dL(1)*r*ss(md(1),md(2),md(3),2)/d2
   END IF
END DO
sys%ss(:,:,:,:) = ss(:,:,:,:)
sys%ss_avg(:,:,:,:) = sys%ss_avg(:,:,:,:) + ss(:,:,:,:)

!
RETURN
END SUBROUTINE HEINZ_STRESS
!
!############################################################################
SUBROUTINE PRINT_HEINZ_stress(sys, Nloop)
USE DATAFMT
!
TYPE(STAT):: sys
INTEGER(GI):: Nloop
!
REAL(DP):: dx(3), box(3), hbox(3)
INTEGER:: i, j, k, n, N_node, N_elem, N_grid(3), id, jd, kd
CHARACTER*7:: name(9), name_avg(9)
REAL(DP), ALLOCATABLE:: ss(:,:), ss_avg(:,:)
INTEGER(GI), ALLOCATABLE:: elem(:,:)

dx(:) = sys%dx(:)
box(:) = sys%box(:)
hbox(:) = box(:)*0.5D0
N_grid(:) = sys%N_grid(:)
name = RESHAPE((/"ss11", "ss22", "ss33", "ss12","ss23","ss31","ssrr", &
     & "sstt", "sspp"/),(/9/))
name_avg = RESHAPE((/"ss11avg", "ss22avg", "ss33avg", "ss12avg","ss23avg", &
     & "ss31avg", "ssrravg", "ssttavg","ssppavg"/),(/9/))
N_node = (N_grid(1) + 1)*(N_grid(2) + 1)*(N_grid(3) + 1)
N_elem = N_grid(1)*N_grid(2)*N_grid(3)
ALLOCATE(ss(N_node,9), ss_avg(N_node,9), elem(N_elem,8), STAT=i)
n = 0
DO k=1, N_grid(3)
   kd = N_grid(2)*N_grid(1)*(k-1)
   DO j=1, N_grid(2)
      jd = N_grid(1)*(j-1)
      DO i=1, N_grid(1)
         n = n + 1
         kd = (N_grid(2)+1)*(N_grid(1)+1)*(k-1)
         jd = (N_grid(1)+1)*(j-1)
         elem(n, 1) = kd + jd + (i-1)
         elem(n, 2) = kd + jd + i
         jd = (N_grid(1)+1)*j
         elem(n, 3) = kd + jd + i
         elem(n, 4) = kd + jd + (i-1)
         kd = (N_grid(2)+1)*(N_grid(1)+1)*k
         jd = (N_grid(1)+1)*(j-1)
         elem(n,5) = kd + jd + (i-1)
         elem(n,6) = kd + jd + i
         jd = (N_grid(1)+1)*j
         elem(n,7) = kd + jd + i
         elem(n,8) = kd + jd + (i-1)
      END DO
   END DO
END DO
!
!
OPEN(UNIT=55, FILE="ss_view.vtk")
!
! VTK file header
WRITE(55,10)
WRITE(55,20)
WRITE(55,30)
10 FORMAT("# vtk DataFile Version 2.0")
20 FORMAT("REAXFF many-body potential stress")
30 FORMAT("ASCII")
!
! Node information
WRITE(55,40)
WRITE(55,50) N_node
n = 0
DO k = 1, N_grid(3) + 1
   DO j = 1, N_grid(2) + 1
      DO i = 1, N_grid(1) + 1
         WRITE(55,60) DBLE(i-1)*dx(1) - hbox(1), DBLE(j-1)*dx(2) - hbox(2), &
              & DBLE(k-1)*dx(3) - hbox(3)
         n = n + 1
         id = i; jd = j; kd = k
         IF (id > N_grid(1)) id = 1
         IF (jd > N_grid(2)) jd = 1
         IF (kd > N_grid(3)) kd = 1
!         PRINT '(7(I3,1X))', N_grid(:), n, id, jd, kd
         ss(n,:) = sys%ss(id,jd,kd,:)
         ss_avg(n,:) = sys%ss_avg(id,jd,kd,:)/DBLE(Nloop)
      END DO
   END DO
END DO
40 FORMAT("DATASET UNSTRUCTURED_GRID")
50 FORMAT("POINTS ", I8,  " float")
60 FORMAT(3(F8.2,1X))
!
! Element structure
WRITE(55, 70) N_elem, 9*N_elem
DO i=1, N_elem
   WRITE(55,80) elem(i,1:8)
END DO
70 FORMAT("CELLS ", I8, 1X, I8)
80 FORMAT("8 ", 8(I8, 1X))
!
! Element type
WRITE(55, 90) N_elem
n = 0
DO i=1, N_elem
   IF (n < 7) THEN   
      WRITE(55,100,ADVANCE="NO")
      n = n + 1
   ELSE 
      WRITE(55,100)
      n = 0
   END IF
END DO
WRITE(55,*)
90 FORMAT("CELL_TYPES", I8)
100 FORMAT("12 ")
!
! Stress output
WRITE(55,110) N_node
DO k=1, 9
   WRITE(55, 120) name(k)
   WRITE(55, 130)
   n = 0
   DO i=1, N_node
      IF (n < 7) THEN
         WRITE(55,140, ADVANCE="NO") ss(i,k)
         n = n + 1
      ELSE
         WRITE(55, 140) ss(i,1)
         n = 0
      END IF
   END DO
   WRITE(55,*)
END DO
!
! Stress_avg output
DO k=1, 9
   WRITE(55, 120) name_avg(k)
   WRITE(55, 130)
   n = 0
   DO i=1, N_node
      IF (n < 7) THEN
         WRITE(55,140, ADVANCE="NO") ss_avg(i,k)
         n = n + 1
      ELSE
         WRITE(55, 140) ss_avg(i,k)
         n = 0
      END IF
   END DO
   WRITE(55,*)
END DO
110 FORMAT("POINT_DATA ", I8)
120 FORMAT("SCALARS ", A7, " float")
130 FORMAT("LOOKUP_TABLE default")
140 FORMAT(ES9.2, 1X)
!
CLOSE(55)
DEALLOCATE(ss, ss_avg, elem, STAT=i)
!
RETURN
END SUBROUTINE PRINT_HEINZ_stress
!
!############################################################################
SUBROUTINE HEINZ_STRESS2(Nset, sys, param, q)
USE DATAFMT
IMPLICIT NONE
!
TYPE(AMNT):: Nset
TYPE(STAT):: sys
TYPE(ATMC):: param(Nset%Natom)
TYPE(PTCL):: q(Nset%Npt_all)
!
INTEGER(GI):: i, j, k, n, id, nd(3), md(3), N_grid(3)
REAL   (DP):: dx(3), dA(3), hbox(3), dv(3), xm, dL(3), r, r2, xx(3), &
     & ss(sys%N_grid(1),sys%N_grid(2),sys%N_grid(3),9), d, d2, sgn(3), tmp

N_grid(:) = sys%N_grid(:)
dx(:) = sys%dx(:)
dA(:) = sys%dA(:)
hbox(:) = sys%box(:)*0.5D0
ss(:,:,:,:) = 0.0D0
dv(:) = dx(:)*dA(:)
!PRINT *, "N grid = ", N_grid(:), dx(:)
DO i=1, Nset%Npt_all
   nd(:) = INT((q(i)%xx(:) + hbox(:))/dx(:)) + 1
   md(:) = nd(:) + 1
   sgn(:) = 1.0D0
   dL(:) = (q(i)%xx(:) + hbox(:)) -  (nd(:) - 1)*dx(:) 
   dL(:) = dL(:)/dx(:)
   !PRINT *, q(i)%xx(:)
   !PRINT *, nd(:)
   !PRINT *, md(:)
   DO j=1,3
      IF (md(j) > N_grid(j)) md(j) = 1
      IF (q(i)%xv(j) < 0.0D0) sgn(j) = -1.D0
   END DO
   id = q(i)%id
   xm = param(id)%am      
   !
   ! kinetic
   !goto 10
   tmp = sgn(1)*xm*q(i)%xv(1)*q(i)%xv(1)/dv(1)
   ss(nd(1),nd(2),nd(3),1) = ss(nd(1),nd(2),nd(3),1) + tmp*(1.D0-dL(1))
   ss(md(1),nd(2),nd(3),1) = ss(md(1),nd(2),nd(3),1) - tmp*dL(1)
   tmp = xm*q(i)%xv(1)*q(i)%xv(2)/dv(1)
   ss(nd(1),nd(2),nd(3),4) = ss(nd(1),nd(2),nd(3),4) + tmp*(1.D0-dL(1))
   ss(md(1),nd(2),nd(3),4) = ss(md(1),nd(2),nd(3),4) - tmp*dL(1)

   tmp = sgn(2)*xm*q(i)%xv(2)*q(i)%xv(2)/dv(2)
   ss(nd(1),nd(2),nd(3),2) = ss(nd(1),nd(2),nd(3),2) + tmp*(1.D0-dL(2))
   ss(nd(1),md(2),nd(3),2) = ss(nd(1),md(2),nd(3),2) - tmp*dL(2)
   tmp = xm*q(i)%xv(2)*q(i)%xv(3)/dv(2)
   ss(nd(1),nd(2),nd(3),5) = ss(nd(1),nd(2),nd(3),5) + tmp*(1.D0-dL(2))
   ss(nd(1),md(2),nd(3),5) = ss(nd(1),md(2),nd(3),5) - tmp*dL(2)
 
   tmp = sgn(3)*xm*q(i)%xv(3)*q(i)%xv(3)/dv(3)
   ss(nd(1),nd(2),nd(3),3) = ss(nd(1),nd(2),nd(3),3) + tmp*(1.D0-dL(3))
   ss(nd(1),nd(2),md(3),3) = ss(nd(1),nd(2),md(3),3) - tmp*dL(3)
   tmp = xm*q(i)%xv(3)*q(i)%xv(1)/dv(3)
   ss(nd(1),nd(2),nd(3),6) = ss(nd(1),nd(2),nd(3),6) + tmp*(1.D0-dL(3))
   ss(nd(1),nd(2),md(3),6) = ss(nd(1),nd(2),md(3),6) - tmp*dL(3)
   !
   ! Force
10 tmp = q(i)%ff(1)/dA(1)
   ss(nd(1),nd(2),nd(3),1) = ss(nd(1),nd(2),nd(3),1) + tmp*(1.D0-dL(1))
   ss(md(1),nd(2),nd(3),1) = ss(md(1),nd(2),nd(3),1) - tmp*dL(1)
   tmp = q(i)%ff(2)/dA(1)
   ss(nd(1),nd(2),nd(3),4) = ss(nd(1),nd(2),nd(3),4) + tmp*(1.D0-dL(1))
   ss(md(1),nd(2),nd(3),4) = ss(md(1),nd(2),nd(3),4) - tmp*dL(1)

   tmp = q(i)%ff(2)/dA(2)
   ss(nd(1),nd(2),nd(3),2) = ss(nd(1),nd(2),nd(3),2) + tmp*(1.D0-dL(2))
   ss(nd(1),md(2),nd(3),2) = ss(nd(1),md(2),nd(3),2) - tmp*dL(2)
   tmp = q(i)%ff(3)/dA(2)
   ss(nd(1),nd(2),nd(3),5) = ss(nd(1),nd(2),nd(3),5) + tmp*(1.D0-dL(2))
   ss(nd(1),md(2),nd(3),5) = ss(nd(1),md(2),nd(3),5) - tmp*dL(2)

   tmp = q(i)%ff(3)/dA(3)
   ss(nd(1),nd(2),nd(3),3) = ss(nd(1),nd(2),nd(3),3) + tmp*(1.D0-dL(3))
   ss(nd(1),nd(2),md(3),3) = ss(nd(1),nd(2),md(3),3) - tmp*dL(3)
   tmp = q(i)%ff(1)/dA(3)
   ss(nd(1),nd(2),nd(3),6) = ss(nd(1),nd(2),nd(3),6) + tmp*(1.D0-dL(3))
   ss(nd(1),nd(2),md(3),6) = ss(nd(1),nd(2),md(3),6) - tmp*dL(3)
END DO
!
DO i=1, N_grid(1)
   xx(1) = dx(1)*DBLE(i-1) - hbox(1)
   DO j=1, N_grid(2)
      xx(2) = dx(2)*DBLE(j-1) - hbox(2)
      DO k=1, N_grid(3)
         xx(3) = dx(3)*DBLE(k-1) - hbox(3)
         dL(:) = hbox(:)
         DO n=1, 3
            IF (xx(n) < 0.) dL(n) = -dL(n)
         END DO
         dL(:) = xx(:) - dL(:)
         r2 = dL(1)*dL(1) + dL(2)*dL(2) + dL(3)*dL(3)
         d2 = dL(1)*dL(1) + dL(2)*dL(2)
         IF (r2 > 0.1 .AND. d2 > 0.1 ) THEN
            ! http://en.wikipedia.org/wiki/List_of_common_coordinate_transformations
            r = DSQRT(r2)
            d = DSQRT(d2)
            dL(:) = dL(:) / r
            ss(i,j,k,7) =  ss(i,j,k,7) + dL(1)*ss(i,j,k,1) + &
                 & dL(2)*ss(i,j,k,2) + dL(3)*ss(i,j,k,3)
            ss(i,j,k,8) =  ss(i,j,k,8) + dL(1)*dL(3)*ss(i,j,k,1)/d + &
                 & dL(2)*dL(3)*ss(i,j,k,2)/d - d*ss(i,j,k,3)/r2
            ss(i,j,k,9) =  ss(i,j,k,9) - dL(2)*r*ss(i,j,k,1)/d2 + &
                 & dL(1)*r*ss(i,j,k,2)/d2
         END IF
      END DO
   END DO
END DO
sys%ss(:,:,:,:) = ss(:,:,:,:)
sys%ss_avg(:,:,:,:) = sys%ss_avg(:,:,:,:) + ss(:,:,:,:)

!
RETURN
END SUBROUTINE HEINZ_STRESS2
!
!############################################################################
SUBROUTINE HEINZ_STRESS3(Nset, sys, param, q)
USE DATAFMT
!
! Using LAMMPS algorithm: http://lammps.sandia.gov/doc/compute_pressure.html
! 
IMPLICIT NONE
!
TYPE(AMNT):: Nset
TYPE(STAT):: sys
TYPE(ATMC):: param(Nset%Natom)
TYPE(PTCL):: q(Nset%Npt_all)
!
INTEGER(GI):: i, j, k, n, id, nd(3), md(3), N_grid(3)
REAL   (DP):: dx(3), dA(3), hbox(3), dv(3), xm, dL(3), r, r2, xx(3), frac(8), &
     & ss(sys%N_grid(1),sys%N_grid(2),sys%N_grid(3),9), d, d2, sgn(3), tmp

N_grid(:) = sys%N_grid(:)
dx(:) = sys%dx(:)
dA(:) = sys%dA(:)
hbox(:) = sys%box(:)*0.5D0
ss(:,:,:,:) = 0.0D0
dv(:) = dx(:)*dA(:)
!PRINT *, "N grid = ", N_grid(:), dx(:)
DO i=1, Nset%Npt_all
   nd(:) = INT((q(i)%xx(:) + hbox(:))/dx(:)) + 1
   md(:) = nd(:) + 1
   sgn(:) = 1.0D0
   dL(:) = (q(i)%xx(:) + hbox(:)) -  (nd(:) - 1)*dx(:) 
   dL(:) = dL(:)/dx(:)
   !PRINT *, q(i)%xx(:)
   !PRINT *, nd(:)
   !PRINT *, md(:)
   DO j=1,3
      IF (md(j) > N_grid(j)) md(j) = 1
      IF (q(i)%xv(j) < 0.0D0) sgn(j) = -1.D0
   END DO
   id = q(i)%id
   xm = param(id)%am      
   frac(1) = dL(1)*dL(2)*dL(3)
   frac(2) = (1.D0 - dL(1))*dL(2)*dL(3)
   frac(3) = dL(1)*(1.D0 - dL(2))*dL(3)
   frac(4) = dL(1)*dL(2)*(1.D0 - dL(3))
   frac(5) = (1.D0 - dL(1))*(1.D0 - dL(2))*dL(3)
   frac(6) = (1.D0 - dL(1))*dL(2)*(1.D0 - dL(3))
   frac(7) = dL(1)*(1.D0 - dL(2))*(1.D0 - dL(3))
   frac(8) = (1.D0 - dL(1))*(1.D0 - dL(2))*(1.D0 - dL(3))
   !
   ! kinetic
   !goto 10
   tmp = q(i)%xv(1)*q(i)%xv(1)/dv(1) + q(i)%xv(2)*q(i)%xv(2)/dv(2) + &
        & q(i)%xv(3)*q(i)%xv(3)/dv(3) + q(i)%xv(1)*q(i)%xv(2)/dv(1) + &
        & q(i)%xv(2)*q(i)%xv(3)/dv(2) + q(i)%xv(3)*q(i)%xv(1)/dv(3)
   tmp = -xm * tmp 
   ss(nd(1),nd(2),nd(3),1) = ss(nd(1),nd(2),nd(3),1) + frac(8)*tmp
   ss(md(1),nd(2),nd(3),1) = ss(md(1),nd(2),nd(3),1) + frac(7)*tmp
   ss(nd(1),md(2),nd(3),1) = ss(nd(1),md(2),nd(3),1) + frac(6)*tmp
   ss(nd(1),nd(2),md(3),1) = ss(nd(1),nd(2),md(3),1) + frac(5)*tmp
   ss(md(1),md(2),nd(3),1) = ss(md(1),md(2),nd(3),1) + frac(4)*tmp
   ss(md(1),nd(2),md(3),1) = ss(md(1),nd(2),md(3),1) + frac(3)*tmp
   ss(nd(1),md(2),md(3),1) = ss(nd(1),md(2),md(3),1) + frac(2)*tmp
   ss(md(1),md(2),md(3),1) = ss(md(1),md(2),md(3),1) + frac(1)*tmp
   !
   ! Force
   ! dA(1) fraction
10 tmp = q(i)%ff(1)/dA(1) + q(i)%ff(2)/dA(1)
   ss(nd(1),nd(2),nd(3),1) = ss(nd(1),nd(2),nd(3),1) + &
        & tmp*(1.D0 - dL(2))*(1.D0 - dL(3))
   ss(md(1),nd(2),nd(3),1) = ss(md(1),nd(2),nd(3),1) - &
        & tmp*(1.D0 - dL(2))*(1.D0 - dL(3))
   ss(nd(1),md(2),nd(3),1) = ss(nd(1),md(2),nd(3),1) + &
        & tmp*dL(2)*(1.D0 - dL(3))
   ss(md(1),md(2),nd(3),1) = ss(md(1),md(2),nd(3),1) - &
        & tmp*dL(2)*(1.D0 - dL(3))
   ss(nd(1),nd(2),md(3),1) = ss(nd(1),nd(2),md(3),1) + &
        & tmp*(1.D0 - dL(2))*dL(3)
   ss(md(1),nd(2),md(3),1) = ss(md(1),nd(2),md(3),1) - &
        & tmp*(1.D0 - dL(2))*dL(3)
   ss(nd(1),md(2),md(3),1) = ss(nd(1),md(2),md(3),1) + &
        & tmp*dL(2)*dL(3)
   ss(md(1),md(2),md(3),1) = ss(md(1),md(2),md(3),1) - &
        & tmp*dL(2)*dL(3)
   !
   ! dA(2) fraction
   tmp = q(i)%ff(2)/dA(2) + q(i)%ff(3)/dA(2)
   ss(nd(1),nd(2),nd(3),1) = ss(nd(1),nd(2),nd(3),1) + &
        & tmp*(1.D0 - dL(1))*(1.D0 - dL(3))
   ss(nd(1),md(2),nd(3),1) = ss(nd(1),md(2),nd(3),1) - &
        & tmp*(1.D0 - dL(1))*(1.D0 - dL(3))
   ss(md(1),nd(2),nd(3),1) = ss(md(1),nd(2),nd(3),1) + &
        & tmp*dL(1)*(1.D0 - dL(3))
   ss(md(1),md(2),nd(3),1) = ss(md(1),md(2),nd(3),1) - &
        & tmp*dL(1)*(1.D0 - dL(3))
   ss(nd(1),nd(2),md(3),1) = ss(nd(1),nd(2),md(3),1) + &
        & tmp*(1.D0 - dL(1))*dL(3)
   ss(nd(1),md(2),md(3),1) = ss(nd(1),md(2),md(3),1) - &
        & tmp*(1.D0 - dL(1))*dL(3)
   ss(md(1),nd(2),md(3),1) = ss(md(1),nd(2),md(3),1) + &
        & tmp*dL(1)*dL(3)
   ss(md(1),md(2),md(3),1) = ss(md(1),md(2),md(3),1) - &
        & tmp*dL(1)*dL(3)
  !
   ! dA(3) fraction
   tmp = q(i)%ff(3)/dA(3) + q(i)%ff(1)/dA(3)
   ss(nd(1),nd(2),nd(3),1) = ss(nd(1),nd(2),nd(3),1) + &
        & tmp*(1.D0 - dL(1))*(1.D0 - dL(2))
   ss(nd(1),nd(2),md(3),1) = ss(nd(1),nd(2),md(3),1) - &
        & tmp*(1.D0 - dL(1))*(1.D0 - dL(2))
   ss(md(1),nd(2),nd(3),1) = ss(md(1),nd(2),nd(3),1) + &
        & tmp*dL(1)*(1.D0 - dL(2))
   ss(md(1),nd(2),md(3),1) = ss(md(1),nd(2),md(3),1) - &
        & tmp*dL(1)*(1.D0 - dL(2))
   ss(nd(1),md(2),nd(3),1) = ss(nd(1),md(2),nd(3),1) + &
        & tmp*(1.D0 - dL(1))*dL(2)
   ss(nd(1),md(2),md(3),1) = ss(nd(1),md(2),md(3),1) - &
        & tmp*(1.D0 - dL(1))*dL(2)
   ss(md(1),md(2),nd(3),1) = ss(md(1),md(2),nd(3),1) + &
        & tmp*dL(1)*dL(3)
   ss(md(1),md(2),md(3),1) = ss(md(1),md(2),md(3),1) - &
        & tmp*dL(1)*dL(3)
END DO
!
DO i=1, N_grid(1)
   xx(1) = dx(1)*DBLE(i-1) - hbox(1)
   DO j=1, N_grid(2)
      xx(2) = dx(2)*DBLE(j-1) - hbox(2)
      DO k=1, N_grid(3)
         xx(3) = dx(3)*DBLE(k-1) - hbox(3)
         dL(:) = hbox(:)
         DO n=1, 3
            IF (xx(n) < 0.) dL(n) = -dL(n)
         END DO
         dL(:) = xx(:) - dL(:)
         r2 = dL(1)*dL(1) + dL(2)*dL(2) + dL(3)*dL(3)
         d2 = dL(1)*dL(1) + dL(2)*dL(2)
         IF (r2 > 0.1 .AND. d2 > 0.1 ) THEN
            ! http://en.wikipedia.org/wiki/List_of_common_coordinate_transformations
            r = DSQRT(r2)
            d = DSQRT(d2)
            dL(:) = dL(:) / r
            ss(i,j,k,7) =  ss(i,j,k,7) + dL(1)*ss(i,j,k,1) + &
                 & dL(2)*ss(i,j,k,2) + dL(3)*ss(i,j,k,3)
            ss(i,j,k,8) =  ss(i,j,k,8) + dL(1)*dL(3)*ss(i,j,k,1)/d + &
                 & dL(2)*dL(3)*ss(i,j,k,2)/d - d*ss(i,j,k,3)/r2
            ss(i,j,k,9) =  ss(i,j,k,9) - dL(2)*r*ss(i,j,k,1)/d2 + &
                 & dL(1)*r*ss(i,j,k,2)/d2
         END IF
      END DO
   END DO
END DO
sys%ss(:,:,:,:) = ss(:,:,:,:)
sys%ss_avg(:,:,:,:) = sys%ss_avg(:,:,:,:) + ss(:,:,:,:)

!
RETURN
END SUBROUTINE HEINZ_STRESS3
