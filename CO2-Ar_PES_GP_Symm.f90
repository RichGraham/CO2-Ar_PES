!!CO2-Ar PES
!!Distances in Angstrom and energies in Hartrees
module PES_CO2_Ar_details
  double precision :: gpRMax = 9.0  
  double precision :: gpRMin = 2.0  
  double precision :: lCO = 1.1632
  double precision :: AngToBohr =  1.8897259885789
  interface PES_CO2_Ar_GP 
     function PES_CO2_Ar_GP(xStar) 
       implicit none 
       double precision:: PES_CO2_Ar_GP
       double precision, dimension(:) ::  xStar
     end function PES_CO2_Ar_GP
  end interface PES_CO2_Ar_GP
end module PES_CO2_Ar_details


module GP_CO2_Ar_variables
  double precision, allocatable :: alpha (:), lScale(:), xTraining(:,:), xTrainingPerm(:,:)
  double precision expVar,NuggVar, gpEmax
  integer :: nDim=3
  integer :: nTraining=221
end module GP_CO2_Ar_variables

double precision function asymp_CO2_Ar(rab)
  use PES_CO2_Ar_details
  implicit none
  double precision rab(3)
  double precision c1,c2,c3, cAr, o1Ar, o2Ar

  c1= (  rab(1)**2 + lCO**2 - rab(2)**2)/2.0/rab(1)/lCO
  
  c2= (  rab(2)**2 + lCO**2 - rab(1)**2)/2.0/rab(2)/lCO
  c3= (  rab(3)**2 + lCO**2 - rab(1)**2)/2.0/rab(3)/lCO


  cAr  = - ( 4.64 * (1 + 3*c1**2)  +  3.30 * (5 - 3*c1**2) ) / (rab(1)*AngToBohr)**6
  o1Ar = - ( 8.69 * (1 + 3*c2**2)  +  4.76 * (5 - 3*c2**2) ) / (rab(2)*AngToBohr)**6
  o2Ar = - ( 8.69 * (1 + 3*c3**2)  +  4.76 * (5 - 3*c3**2) ) / (rab(3)*AngToBohr)**6
  
  asymp_CO2_Ar= cAr + o1Ar + o2Ar
end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Gaussian Process Code!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!pom


  
subroutine load_GP_CO2_Ar_Data
  use GP_CO2_Ar_variables
  use PES_CO2_Ar_details
  implicit none
  
  double precision, allocatable::  xStar(:)
  integer i,j
  double precision :: dum, expVar1, expVar2
  character (len=90) :: filename
  CHARACTER(len=255) :: homedir,codedir

  allocate (alpha(nTraining), lScale(nDim), xTraining(nDim,nTraining),xTrainingPerm(nDim,nTraining), xStar(nDim))
  CALL getenv("HOME", homedir)
  codedir=TRIM(homedir) // '/source/CO2-Ar_PES'
  !call chdir(codedir)
  
  !====Load hyperparameters====
  write (filename, '( "/TrainingData/HyperParams_Symm", I3.3, ".dat" )' )  nTraining
  filename = trim(codedir) // filename
  open (unit = 7, file = filename)
  !Only need to read some as others are tied.
  read (7,*) lScale(2), lScale(1), expVar1, expVar2,NuggVar, gpEmax
  !Copy over the tied values
  lScale(3) = lScale(2)
  expVar = expVar1 * expVar2
  !print *,"HyperParams",lScale(1), lScale(2), lScale(3), expVar1, expVar2,NuggVar, gpEmax
  close(7)
  
  !====Load alpha coefficients====
  write (filename, '( "/TrainingData/alpha_Symm", I3.3, ".dat" )' )  nTraining
  filename = trim(codedir) // filename
  open (unit = 7, file = filename)
  do i=1,nTraining
     read (7,*) alpha(i)
     !print *,"alpha ",i, alpha(i)
  end do
  close(7)


  !====Load training data x values ====
  write (filename, '( "/TrainingData/xTraining", I3.3, ".dat" )' )  nTraining
  filename = trim(codedir) // filename
  open (unit = 7, file = filename)
    
  do i=1,nTraining
     read (7,*) xTraining(1,i), xTraining(2,i), xTraining(3,i)
     !print *,xTraining(1,i), xTraining(2,i), xTraining(3,i)
  end do
  !close(7)
  
  !! Permute the training vectors
  xTrainingPerm = xTraining
  do i=1,nTraining
     xTrainingPerm(2,i)=xTraining(3,i)
     xTrainingPerm(3,i)=xTraining(2,i)
  end do

end subroutine load_GP_CO2_Ar_Data
  
function PES_CO2_Ar_GP(xStar)
  use GP_CO2_Ar_variables
  implicit none
  double precision, dimension(:) :: xStar
  double precision:: PES_CO2_Ar_GP
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
  
  PES_CO2_Ar_GP=kKern * expVar
end function PES_CO2_Ar_GP



function PES_CO2_Ar( rab )
  !! Takes in rab in Angstrom
  use PES_CO2_Ar_details
  use GP_CO2_Ar_variables
  implicit none
  double precision rab(3), xStar(3), asymp_CO2_Ar
  double precision  PES_CO2_Ar
  double precision repFactor

  repFactor=1.0
  
  if( rab(1) > gpRMax  .AND.  rab(2) > gpRMax .AND.  rab(3) > gpRMax &
       ) then !!Use asymptotic function
     PES_CO2_Ar = asymp_CO2_Ar(rab)
     
  else if (rab(1) < gpRMin/repFactor  .OR.  rab(2) < gpRMin/repFactor  .OR.  rab(3) < gpRMin/repFactor &
       ) then !! Use repulsive approximation function
     PES_CO2_Ar=gpEmax* (1.0/rab(1)**12+1.0/rab(2)**12+1.0/rab(3)**12) *gpRMin **12
     
  else !! Use the Guassian Process function
     xStar(:) = 1/rab(:)
     PES_CO2_Ar = PES_CO2_Ar_GP( xStar)
  end if

  !PES=gpEmax/3.0* (1.0/rab(1)**12+1.0/rab(2)**12+1.0/rab(3)**12) *gpRMin **12

  
end function PES_CO2_Ar
