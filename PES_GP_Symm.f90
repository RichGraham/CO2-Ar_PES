module PES_details
  double precision :: gpRMax = 8.0   *   1.8897259885789
  double precision :: gpRMin = 2.0   *   1.8897259885789
  double precision :: gpEmax = 0.456
  !! Maximum energy in Hartree seen within the geometric constraint
  interface PES_GP 
     function PES_GP(xStar) 
       implicit none 
       double precision:: PES_GP
       double precision, dimension(:) ::  xStar
     end function PES_GP
  end interface PES_GP
end module PES_details


module GP_variables
  double precision, allocatable :: alpha (:), lScale(:), xTraining(:,:), xTrainingPerm(:,:)
  double precision expVar,NuggVar
  integer :: nDim=6
  integer :: nTraining=135 
end module GP_variables


 ! Test program
implicit none
double precision rco2, rco
double precision angDat(5)

integer k,i, choice

! Bond length of CO, NB, Bohr
rco=2.132d0

! Bond length of CO2, NB, Bohr
rco2=2.198d0

call load_GP_Data
call fixedAngleSlice(rco2,rco)

end
!

subroutine fixedAngleSlice( rco2, rco)
  use PES_details
  implicit none
  double precision rco2, rco,rab(6)
  integer i, itot
  double precision  r, tha, thb, ph, e, e_GP, asymp, PES,AngToBohr
  AngToBohr= 1.8897259885789
    
  itot=500
  tha =  3.14159265359/2.0
  thb =  3.14159265359
  ph =  0.0

  open (unit=15, file="PES_Out.dat ", status='replace')
  
  do i=0, itot

     ! specify centre-to-centre separation in Bohr
     r = (  2.4 + 9.0*i/(1.0*itot) ) * AngToBohr

     call computeDistances(r,tha,thb,ph,rab, rco2, rco)
     
     
     e=PES( rab, rco2, rco)
     !e_GP = PES_GP( xStar)
     write(15,*) r/AngToBohr , e 
     
  enddo
  write(6,*)'Written to file: PES_Out.dat '
  close(15)

end subroutine fixedAngleSlice
  
  
  
subroutine computeDistances(r,tha,thb, ph, rab,rco2,rco)
  implicit none
  double precision rco2, rco, rab(6), ca(3,3), cb(3,2), axisa(3), axisb(3)
  double precision centa(3), centb(3), r, tha, thb, ph,   AngToBohr
  integer ia, ib, ir,k

  AngToBohr= 1.8897259885789
     
!!! Load positions of atoms
   !====Position CO2====
   ! Direction and centre of CO2
   !1=x; 2=y; 3=y
   centa(:)=0d0
   axisa(1)=sin(tha)
   axisa(2)=0d0
   axisa(3)=cos(tha)

   ! Atom positions of CO2
   ca(:,1)=centa(:)+0d0 ! C1
   ca(:,2)=centa(:)+axisa(:)*rco2 ! O1
   ca(:,3)=centa(:)-axisa(:)*rco2 ! O2
   
   
   !====Position CO====
   ! Direction and centre of CO
   centb(:)=0d0
   centb(3)=r
   
   axisb(1)=sin(thb)*cos(ph)
   axisb(2)=sin(thb)*sin(ph)
   axisb(3)=cos(thb)

   ! Atom positions of CO, each half a bond length from centre.
   
   !====RSG Change====
   !Changed the definition of axis so that O is along the minus axis and C is along +ve axis
   !Needed to be consistent with mol pro angles - does not affect final calculation as distances (not angles) are used
   !cb(:,1)=centb(:)+rco*axisb(:)/2d0 ! O'
   !cb(:,2)=centb(:)-rco*axisb(:)/2d0 ! C'
   cb(:,1)=centb(:)-rco*axisb(:)/2d0 ! O'
   cb(:,2)=centb(:)+rco*axisb(:)/2d0 ! C'
   

   
   !Compute interatomic distances
   ! Work out the interatomic distances CO', CC', OO', OC', O-O', O-C'
   ir=1
   do ia=1,3 !atom i of molecule a
      do ib=1,2 !atom i of molecule b
         
         !!Compute the distance between atom ia and ib
         rab(ir)=0
         do k=1,3 !! Cartesian component
            rab(ir)=rab(ir)+(cb(k,ib)-ca(k,ia))**2
         end do
         rab(ir)=sqrt(rab(ir))
         ir=ir+1
      end do
   end do

   
   
end subroutine

double precision function asymp(rab,rco2,rco)
! Work out the asymptotic energy for CO2-CO
! rab are the interatomic distances CO', CC', OO', OC', O-O', O-C' <-These are consistent with RSG vector of distances r1....r6 used in Python and mol pro conversion/lhc codes
! rco2 and rco are CO distances in the molecules
! Everything is in Atomic Units (Bohr, Hartree)
implicit none
double precision rab(6),rco2,rco,emult
double precision qc,qo,qcp,qop,uc,uo,ucp,uop,tcp,top
double precision cozz,cczz,cozx,cczx,coxz,ccxz,coxx,ccxx
double precision oozz,oczz,oozx,oczx,ooxz,ocxz,ooxx,ocxx
! The energy consists of six pair contributions, each containing
! seven terms (six electrostatic and one dispersion).
! The function emult is used to get the seven terms.
! This function calls emult six times, and adds the results.
asymp=0d0
! The arguments of emult are the relevant coordinates and parameters.
! Coordinates are rab, rab', ra'b, ra'b', raa', rbb' where rab is the distance
! between the atoms of interest, a' is bonded to a, b' is bonded to b.
! Parameters are charge of a, dipole of a (in the direction away from
! a'), charge of b, dipole of b (away from b'), quadrupole of b,
! and four ab dispersion energy coefficients, ||, |-, -| and --
! where | means parallel to the bond and - means perpendicular.
!
! Parameters
qc=0.8670d0 ! charges
qo=-0.4335d0
qcp=0.061d0
qop=-0.061d0
uc=0d0 ! dipoles; C doesn't have one
uo=0.1347d0 ! away from C
ucp=-0.186d0 ! away from O'
uop=-0.008d0 ! away from C'
tcp=-0.743d0 ! quadrupoles, CO only
top=-0.189d0 ! quadrupoles, CO only
cozz=2.66d0 ! dispersion, C of CO2, O of CO, z means parallel
cczz=3.45d0
cozx=1.42d0 ! x means perpendicular
cczx=3.26d0
coxz=1.88d0
ccxz=2.44d0
coxx=1.01d0
ccxx=2.31d0
oozz=4.96d0
oczz=6.44d0
oozx=2.66d0
oczx=6.11d0
ooxz=2.70d0
ocxz=3.51d0
ooxx=1.45d0
ocxx=3.34d0
!
! CO' (so a=C, b=O', a'=O-, b'=C')
asymp=asymp+emult(rab(1),rab(2),rab(5),rab(6),rco2,rco,qc,uc,qop,uop,top,cozz,cozx,coxz,coxx)
! CC' (so a=C, b=C', a'=O-, b'=O')
asymp=asymp+emult(rab(2),rab(1),rab(6),rab(5),rco2,rco,qc,uc,qcp,ucp,tcp,cczz,cczx,ccxz,ccxx)
! OO' (so a=O, b=O', a'=C, b'=C')
asymp=asymp+emult(rab(3),rab(4),rab(1),rab(2),rco2,rco,qo,uo,qop,uop,top,oozz,oozx,ooxz,ooxx)
! OC' (so a=O, b=C', a'=C, b'=O')
asymp=asymp+emult(rab(4),rab(3),rab(2),rab(1),rco2,rco,qo,uo,qcp,ucp,tcp,oczz,oczx,ocxz,ocxx)
! O-O' (so a=O-, b=O', a'=C, b'=C')
asymp=asymp+emult(rab(5),rab(6),rab(1),rab(2),rco2,rco,qo,uo,qop,uop,top,oozz,oozx,ooxz,ooxx)
! O-C' (so a=O-, b=C', a'=C, b'=O')
asymp=asymp+emult(rab(6),rab(5),rab(2),rab(1),rco2,rco,qo,uo,qcp,ucp,tcp,oczz,oczx,ocxz,ocxx)
end
!
double precision function emult(rab,rabp,rapb,rapbp,raap,rbbp,qa,ua,qb,ub,tb,d1,d2,d3,d4)
implicit none
double precision rab,rabp,rapb,rapbp,raap,rbbp,qa,ua,qb,ub,tb,d1,d2,d3,d4
double precision eadd, costa, costb, cosphi, dc1, dc2, dc3, dc4
emult=0d0
! Calculate the multipolar energy (Coulomb+dispersion)
!write(6,*)'emult called with distances ',rab,rabp,rapb,rapbp,raap,rbbp
!write(6,*)'and parameters ',qa,ua,qb,ub,tb,d1,d2,d3,d4
! Charges
eadd=qa*qb/rab
!write(6,*)'Charge-charge ',eadd
emult=emult+eadd
! Dipole of A with charge of B, for which we need the angle ta between
! the ap->a axis and the a->b vector
costa=(rapb**2-raap**2-rab**2)/(2d0*raap*rab)
!write(6,*)'cos(ta)=',costa
eadd=ua*qb*costa/rab**2
!write(6,*)'Dipole-charge ',eadd
emult=emult+eadd
! Charge of A with dipole of B, for which we need the angle tb between
! the bp->b axis and the b->a vector
costb=(rabp**2-rbbp**2-rab**2)/(2d0*rbbp*rab)
!write(6,*)'cos(tb)=',costb
eadd=qa*ub*costb/rab**2
!write(6,*)'Charge-dipole ',eadd
emult=emult+eadd
! Dipole of A with dipole of B, for which we need the angle phi between
! the ap->a axis and the bp->b axis (NB, not usual definition of phi)
cosphi=(rabp**2+rapb**2-rab**2-rapbp**2)/(2d0*raap*rbbp)
!write(6,*)'cos(phi)=',cosphi
eadd=ua*ub*(cosphi+3d0*costa*costb)/rab**3
!write(6,*)'Dipole-dipole ',eadd
emult=emult+eadd
! Charge of A with quadrupole of B
eadd=qa*tb*(3d0*costb*costb-1d0)/(2d0*rab**3)
!write(6,*)'Charge-quadrupole ',eadd
emult=emult+eadd
! Dipole of A with quadrupole of B
eadd=ua*tb*(15d0*costa*costb*costb-3d0*costa+6d0*costb*cosphi)/(2d0*rab**4)
!write(6,*)'Dipole-quadrupole ',eadd
emult=emult+eadd
! Dispersion par-par
dc1=d1*(cosphi+3*costa*costb)**2
! Dispersion par-perp
dc2=d2*((1+3*costa**2)-(cosphi+3*costa*costb)**2)
! Dispersion perp-par
dc3=d3*((1+3*costb**2)-(cosphi+3*costa*costb)**2)
! Dispersion perp-perp
dc4=d4*((4-3*costa**2-3*costb**2)+(cosphi+3*costa*costb)**2)
!write(6,*)'DCs ',dc1,dc2,dc3,dc4
eadd=-(dc1+dc2+dc3+dc4)/rab**6
!write(6,*)'Dispersion ',eadd
emult=emult+eadd
!write(6,*)'Total ',emult
end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Gaussian Process Code!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!pom


  
subroutine load_GP_Data
  use GP_variables
  use PES_details
  implicit none
  
  double precision, allocatable::  xStar(:)
  integer i,j
  double precision :: dum
  character (len=90) :: filename

  allocate (alpha(nTraining), lScale(nDim), xTraining(nDim,nTraining),xTrainingPerm(nDim,nTraining), xStar(nDim))

  !====Load hyperparameters====
  write (filename, '( "N", I4.4, "/HyperParams_Symm.dat" )' )  nTraining
  open (unit = 7, file = filename)
  !Only need to read some as others are tied.
  read (7,*) lScale(4), lScale(3), lScale(2), lScale(1),expVar,NuggVar
  !Copy over the tied values
  lScale(5) = lScale(3)
  lScale(6) = lScale(4)
  !print *,"HyperParams",lScale(1), lScale(2), lScale(3), lScale(4), lScale(5), lScale(6),expVar,NuggVar
  close(7)
  
  !====Load alpha coefficients====
  write (filename, '( "N", I4.4, "/alpha_Symm.dat" )' )  nTraining
  open (unit = 7, file = filename)
  do i=1,nTraining
     read (7,*) alpha(i)
     !print *,"alpha ",i, alpha(i)
  end do
  close(7)


  !====Load training data x values ====
  write (filename, '( "N", I4.4, "/xTraining.dat" )' )  nTraining
  open (unit = 7, file = filename)
    
  do i=1,nTraining
     read (7,*) xTraining(1,i), xTraining(2,i), xTraining(3,i), xTraining(4,i), xTraining(5,i), xTraining(6,i)  
  end do
  close(7)
  
  !! Permute the training vectors
  xTrainingPerm = xTraining
  do i=1,nTraining
     xTrainingPerm(3,i)=xTraining(5,i)
     xTrainingPerm(5,i)=xTraining(3,i)
     xTrainingPerm(4,i)=xTraining(6,i)
     xTrainingPerm(6,i)=xTraining(4,i)
  end do

end subroutine load_GP_Data
  
function PES_GP(xStar)
  use GP_variables
  implicit none
  double precision, dimension(:) :: xStar
  double precision:: PES_GP
  integer beta,i
  double precision kSqExp, kSqExpPerm, kKern

  kKern=0

  do i=1,nTraining
     kSqExp=1;
     kSqExpPerm=1;
     do beta=1,nDim
        kSqExp  =  kSqExp * ( exp( - (xStar(beta)-xTraining(beta,i))**2 /2.0/lScale(beta)**2) )
        kSqExpPerm  =  kSqExpPerm * ( exp( - (xStar(beta)-xTrainingPerm(beta,i))**2 /2.0/lScale(beta)**2) ) 
     end do
     kKern = kKern + alpha(i) * (kSqExp + kSqExpPerm)
  end do
  
  PES_GP=kKern * expVar
end function PES_GP



function PES( rab, rco2, rco )
  !! Takes in rab in Angstrom
  use PES_details
  implicit none
  double precision rab(6), xStar(6), rco2, rco, asymp
  double precision AngToBohr, PES

  AngToBohr= 1.8897259885789
  
  if( rab(1) > gpRMax  .AND.  rab(2) > gpRMax .AND.  rab(3) > gpRMax .AND. &
       rab(4) > gpRMax .AND.  rab(5) > gpRMax .AND.  rab(5) > gpRMax  &
       ) then !!Use asymptotic function
     PES = asymp(rab,rco2,rco)
     
  else if (rab(1) < gpRMin  .OR.  rab(2) < gpRMin  .OR.  rab(3) < gpRMin  .OR.  &
       rab(4) < gpRMin  .OR.  rab(5) < gpRMin  .OR.  rab(6) < gpRMin  &
       ) then !! Use repulsive approximation function
     PES=gpEmax/6.0* ( 1.0/rab(1)**12 +   1.0/rab(2)**12 +  1.0/rab(3)**12       &
          + 1.0/rab(4)**12  +   1.0/rab(5)**12  +   1.0/rab(6)**12 ) * gpRMin **12
     
  else !! Use the Guassian Process function
     xStar(:) = AngToBohr/rab(:)
     PES = PES_GP( xStar)
  end if

end function PES
