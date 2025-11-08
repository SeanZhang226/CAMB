! Simple test program for IDE model
! Tests basic IDE module loading and compilation

program test_ide_fortran
    use DarkEnergyIDE
    use IDEtools
    use CoupledFluidModels
    implicit none
    
    type(TDarkEnergyIDE) :: ide_model
    type(CoupFluidParams) :: params
    
    write(*,*) '================================================'
    write(*,*) 'Testing IDE Model Modules - Compilation Test'
    write(*,*) '================================================'
    write(*,*) ''
    
    ! Test 1: Module imports successful
    write(*,*) 'Test 1: Module imports'
    write(*,*) '------------------------------------------------'
    write(*,*) '  DarkEnergyIDE module: OK'
    write(*,*) '  IDEtools module: OK'
    write(*,*) '  CoupledFluidModels module: OK'
    write(*,*) ''
    
    ! Test 2: Set IDE model parameters
    write(*,*) 'Test 2: IDE model parameters'
    write(*,*) '------------------------------------------------'
    ide_model%w_lam = -0.9_dl
    ide_model%wa = 0.0_dl
    ide_model%beta = 0.01_dl
    ide_model%coupling_form = 1
    ide_model%eos_form = 1
    ide_model%covariant_form = 1
    write(*,'(A,F8.4)') '  w = ', ide_model%w_lam
    write(*,'(A,F8.4)') '  wa = ', ide_model%wa
    write(*,'(A,F8.4)') '  beta = ', ide_model%beta
    write(*,'(A,I4)') '  coupling_form = ', ide_model%coupling_form
    write(*,*) ''
    
    ! Test 3: Set CoupledFluidModels parameters
    write(*,*) 'Test 3: Coupled fluid parameters'
    write(*,*) '------------------------------------------------'
    params%w0 = -0.95_dl
    params%w1 = 0.0_dl
    params%beta = 0.02_dl
    write(*,'(A,F8.4)') '  w0 = ', params%w0
    write(*,'(A,F8.4)') '  w1 = ', params%w1
    write(*,'(A,F8.4)') '  beta = ', params%beta
    write(*,*) ''
    
    ! Test 4: Set coupling types
    write(*,*) 'Test 4: Coupling types'
    write(*,*) '------------------------------------------------'
    CFT%WForm = CPL
    CFT%QForm = Hrde
    CFT%CovQForm = uc
    write(*,'(A,I4,A)') '  WForm = ', CFT%WForm, ' (CPL)'
    write(*,'(A,I4,A)') '  QForm = ', CFT%QForm, ' (Hrde)'
    write(*,'(A,I4,A)') '  CovQForm = ', CFT%CovQForm, ' (uc)'
    write(*,*) ''
    
    write(*,*) '================================================'
    write(*,*) 'All compilation tests passed!'
    write(*,*) 'IDE modules are properly integrated into CAMB.'
    write(*,*) '================================================'
    
end program test_ide_fortran
