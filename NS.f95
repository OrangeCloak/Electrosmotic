program NS
    implicit none
    
    ! -------------------------------------------- Declaring Variables : ------------------------------------------------------------
    ! Defining Domain 
    integer :: nx , ny ,ir , jr
    double precision :: dx , dy ,lx , ly
    
    ! Define Flow Properties
    double precision :: Re , G , alpha , beta
    
    ! Define Time Constraints
    real :: dt
    integer :: nt

    ! Define Coordinates Array


    ! Define Velocity and Pressure Array 
    double precision :: uold(101,21) , ustar(101,21) , unew(101,21) ,ufinal(101,21)
    double precision :: vold(101,21) , vstar(101,21) , vnew(101,21) ,vfinal(101,21)
    double precision :: pold(101,21) , pstar(101,21) , pnew(101,21) , pcor(101,21) , bp(101,21) ,pfinal(101,21)
    double precision :: phi(101,21) , phif(101,21)


    ! Define Iteration Parameters
    integer :: iter
    integer :: time_iter
    double precision :: error , verror


    ! Define SOR Parameters 
    double precision :: n1 , n2 , n3 , n4
    double precision , dimension(101,21) :: ue , uw  
    double precision , dimension(101,21) :: vs , vn 
    double precision , dimension(101,21) :: ape , apw , aps , apn , app
    double precision , dimension(101,21) :: ae , aw , an , as , ap
    double precision , dimension(101,21) :: phiE , phiW , phiN , phiS ,phiTotal , A, B
    integer :: i ,j 



    !-------------------------------------------------------------------------------------------------------------------------------



    ! Declare Domain Constraints 

    nx = 101
    ny = 21
    ir = nx - 1
    jr = ny - 1
    lx = 5
    ly = 1

    dx = lx/ir
    dy = ly/jr

    Re = 0.0427
    G = 23.42831211
    alpha = 2.435560743
    beta = 3.7e-23


    dt = 0.001
    nt = 10


    ! Initialize Velocities and Pressure Array
    unew(:,:) = 0.0
    uold(:,:) = 0.0
    ustar(:,:) = 0.0
    phi(:,:) = 0.0

    vold(:,:) = 0.0
    vstar(:,:) = 0.0
    vnew(:,:) = 0.0

    pstar(:,:) = 1.0d0
    pold(:,:) = 0.0
    pnew(:,:) = 0.0
    pcor(:,:) = 0.0
    bp(:,:) = 0.0
    
    phi(:,:) = 0.0

    ! Inititialize Error 
    error = 1.0d0
    iter = 0
    time_iter = 1






    !---------------------------------------------            S.O.R            ----------------------------------------------------

    do while (time_iter < nt )
        do while (error > 1.0E-06)


            ! X - Momentum Equation
            n1=(dt)*(1.0d0/(2*dx))
            n2=(dt)*(1.0d0/(2*dy))
            n3=(1.0d0/(dx**2))*dt*(1.0d0/Re)
            n4=(1.0d0/(dy**2))*dt*(1.0d0/Re)

            do i = 2, (ir-1)
                do j = 2, jr

                    ue(i,j)=0.5*(uold(i,j)+uold(i+1,j))
                    ae(i,j)=n3 - n1*ue(i,j)
                    uw(i,j)=0.5*(uold(i,j)+uold(i-1,j))
                    aw(i,j)=n3 + n1*uw(i,j)
                    vn(i,j)=0.5*(vold(i,j)+vold(i+1,j))
                    an(i,j)=n4 - n2*vn(i,j)
                    vs(i,j)=0.5*(vold(i,j-1)+vold(i+1,j-1))
                    as(i,j)=n4 + n2*vs(i,j)

                    ap(i,j)= 1 + n1*(uw(i,j)-ue(i,j))+ n2*(vs(i,j)-vn(i,j))- 2*(n3+n4)

                    ! Define Potential Term :
                    phiE(i,j) = 0.5*(phi(i+1,j)-phi(i,j))
                    phiW(i,j) = 0.5*(phi(i,j)-phi(i-1,j))
                    phiN(i,j) = 0.5*(phi(i,j+1)-phi(i,j))
                    phiS(i,j) = 0.5*(phi(i,j)-phi(i,j-1))

                    phiTotal(i,j) = (phiE(i,j) - phiW(i,j))/dx + (phiN(i,j) - phiS(i,j))/dy

                    ustar(i,j)=(ae(i,j)*uold(i+1,j)+aw(i,j)*uold(i-1,j) + an(i,j)*uold(i,j+1)+as(i,j)*uold(i,j-1)+ ap(i,j)*uold(i,j)) + dt*(1.0d0/dx)*(pold(i,j)-pold(i+1,j) + G*phiTotal(i,j)) 
                    
                end do
            end do



            ! X - Equation Boundary Condition
            !Right wall
            do j=1,(ny-1)
                i=(nx-1)
                ustar(i,j)=ustar(i-1,j)
            end do
            
            !Left Wall
            do j=1,ny
                i=1
                ustar(i,j)=1.0d0

            end do
            
            !Top wall
            do i=2,(nx-1)
                j=ny
                ustar(i,j)= -ustar(i,j-1)

            end do
            
            !Bottom Wall
            do i=2,(nx-1)
                j=1
                ustar(i,j)= -ustar(i,j+1)

            end do
                

            ! Y Momentum Equation 
            do i = 2, (nx-1)
                do j = 2, (ny-2)
                    
                    ue(i,j)=0.5*(uold(i,j)+uold(i,j+1))
                    ae(i,j)=n3 - n1*ue(i,j)
                    uw(i,j)=0.5*(uold(i-1,j)+uold(i-1,j+1))
                    aw(i,j)=n3 + n1*uw(i,j)
                    vn(i,j)=0.5*(vold(i,j)+vold(i,j+1))
                    an(i,j)=n4 - n2*vn(i,j)
                    vs(i,j)=0.5*(vold(i,j)+vold(i,j-1))
                    as(i,j)=n4 + n2*vs(i,j)

                    ap(i,j)= 1 + n1*(uw(i,j)-ue(i,j)) + n2*(vs(i,j)-vn(i,j))- 2*(n3+n4)

                    vstar(i,j) = (ae(i,j)*vold(i+1,j)+aw(i,j)*vold(i-1,j) + an(i,j)*vold(i,j+1)+as(i,j)*vold(i,j-1) + ap(i,j)*vold(i,j))+(dt/dy)*(pold(i,j)-pold(i,j+1))

                end do
            end do


            ! Y - Momentum Boundary Condition
            do j=2,(ny-1)
                i=nx
                vstar(i,j)=vstar(i-1,j)
            end do
                    
            do j=2,(ny-1)
                i=1
                vstar(i,j)=0.0
            end do
                
            do i=1,nx
                j=(ny-1)
                vstar(i,j)=0.0
            end do
                
            do i=1,nx
                j=1
                vstar(i,j)=0.0
            end do
            

            !------------------------------                Pressure Correction        ---------------------------------------------------
            pcor(:,:) = 0.0
            do i = 2, (nx-1)
                do j = 2, (ny-1)
                    
                    ape(i,j)=(1.0d0/dy)*dt*dx
                    apw(i,j)=(1.0d0/dy)*dt*dx
                    apn(i,j)=(1.0d0/dx)*dt*dy
                    aps(i,j)=(1.0d0/dx)*dt*dy
                    app(i,j)=ape(i,j)+apw(i,j)+apn(i,j)+aps(i,j)

                    bp(i,j)=((ustar(i-1,j)-ustar(i,j))*dy + (vstar(i,j-1) - vstar(i,j))*dx)

                    pcor(i,j) = (ape(i,j)*pcor(i+1,j) + apw(i,j)*pcor(i-1,j) + apn(i,j)*pcor(i,j+1) + aps(i,j)*pcor(i,j-1)+bp(i,j))/app(i,j)

                end do
            end do
 
            do i=2,(nx-1)
                do j=2,(ny-1)
                    pnew(i,j)=pold(i,j) + 0.8*pcor(i,j)
                end do
            end do
            

            ! Pressure Boundary Condition
            !Top Wall
            do i=1,nx
                j=ny
                pnew(i,j)=pnew(i,j-1)
            end do
            
            !Bottom Wall
            do i=1,nx
                j=1
                pnew(i,j)=pnew(i,j+1)
            end do
            
            !Left Wall
            do j=1,(ny-1)
                i=1
                pnew(i,j)=pnew(i+1,j)
            end do

            !Right Wall
            do j = 1, (ny-1)
                i=nx
                pnew(i,j)=pnew(i-1,j)
            end do





            ! ---------------------------------------       Electric Fields Equation        -------------------------------------



            ! Poisson - Boltzmann Equation :

            phi(:,:) = 0.0
            verror = 0.0

            do while (verror > 1e-06)
                
                do i = 2, (nx-1)
                    do j = 2, (ny-1)
                        
                        A(i,j) = beta*alpha*cosh(alpha*phi(i,j))
                        B(i,j) = beta*sinh(alpha*phi(i,j)) - beta*alpha*phi(i,j)*cosh(alpha*phi(i,j))
    
                        phiE(i,j) = 0.5*(phi(i+1,j)-phi(i,j))
                        phiW(i,j) = 0.5*(phi(i,j)-phi(i-1,j))
                        phiN(i,j) = 0.5*(phi(i,j+1)-phi(i,j))
                        phiS(i,j) = 0.5*(phi(i,j)-phi(i,j-1))
    
                        phi(i,j) = (  ((phiE(i,j) - phiW(i,j))/dy) + ((phiN(i,j) - phiS(i,j))/dx) - B(i,j)  )/A(i,j)
    
                        verror=verror+dabs(phi(i,j))
    
                    end do
                end do

            end do
                

            ! Phi Boundary Conditions

            !Top Wall
            do i=1,nx
                j=ny
                phi(i,j)= 2.0d0 - phi(i,j-1)
            end do
            
            !Bottom Wall
            do i=1,nx
                j=1
                phi(i,j)= 2.0d0 - phi(i,j+1)
            end do
            
            !Left Wall
            do j=1,(ny-1)
                i=1
                phi(i,j)=phi(i+1,j)
            end do

            !Right Wall
            do j = 1, (ny-1)
                i=nx
                phi(i,j)=phi(i-1,j)
            end do





            ! ---------------------------------------        Velocity Correction           ---------------------------------------
            do i=2,(nx-1)
                do j=2,(ny-1)
                    unew(i,j)= ustar(i,j)+dt*(1.0d0/dx)*(pcor(i,j)-pcor(i+1,j))
                end do
            end do

            do i=2,(nx-1)
                do j=2,(ny-1)
                    vnew(i,j)= vstar(i,j)+(dt)*(1.0d0/dy)*(pcor(i,j)-pcor(i,j+1))
                end do
            end do

            ! X - Velocity Boundary Conditions
            !Right Wall
            do j=1,(ny-1)
                i=(nx-1)
                unew(i,j)=unew(i-1,j)
            end do
            
            !Left Wall
            do j=1,ny
                i=1
                unew(i,j)=1.0d0
            end do
            
            !Top Wall
            do i=2,(nx-1)
                j=ny
                unew(i,j)= -unew(i,j-1)
            end do
            
            !Bottom Wall
            do i=2,(nx-1)
                j=1
                unew(i,j)= -unew(i,j+1)
            end do

            ! Y - Velocity Boundary Conditions
            do j=2,(ny-1)
                i=nx
                vnew(i,j)=vnew(i-1,j)
            end do
                    
            do j=2,(ny-1)
                i=1
                vnew(i,j)=0.0
            end do
                
            do i=1,nx
                j=(ny-1)
                vnew(i,j)=0.0
            end do
                
            do i=1,nx
                j=1
                vnew(i,j)=0.0 
            end do


            ! ----------------------------------------             Error Handling            --------------------------------------
            error = 0.0
            do i=1,nx
                do j=1,ny
                    error=error+dabs(bp(i,j))
                end do
            end do

            iter = iter + 1
            print *, "Current error : " ,error , "Current iteration :",iter

            ! Update Velocities and Pressure 
            uold = unew
            vold = vnew
            pold = pnew

        end do 
        error = 1
        time_iter = time_iter + 1
        dt = 0.001 + dt

        print *, "--------------------"
        print *, "Current time loop :" , time_iter , " Time Step : ", dt
        print *, "--------------------"

    end do 

    print *, "Domain spacing is " , dx ,dy
    print *, "Reynold's Number is" , Re
    print *, "Time Step Length :" , dt , "Number Of Time Steps :" ,nt 

    ! -------------------------------------------------      Reallocating Arrays        -------------------------------------------

    ufinal(:,:) = 0.0
    vfinal(:,:) = 0.0
    pfinal(:,:) = 0.0
    phif(:,:) = 0.0

    do i=1,(nx-1)
        do j=1,(ny-1)
            ufinal(i,j)=0.5*(unew(i,j) + unew(i,j+1))
            vfinal(i,j)=0.5*(vnew(i,j) + vnew(i+1,j))
            pfinal(i,j)=0.25*(pnew(i,j) + pnew(i,j+1) + pnew(i+1,j) + pnew(i+1,j+1))
            phif(i,j) = 0.25*(phif(i,j) + phif(i,j+1) + phif(i+1,j) + phif(i+1,j+1))

        end do
    end do
        

    ! Exporting Data to the file :
    open(unit = 2 , file='uniform.csv')
    do i = 1, (nx-1)
        do j = 1, (ny-1)
            write(2,*) ufinal(i,j) , ',' , vfinal(i,j) , ',' , pfinal(i,j) 
        end do
    end do
    close(2)

end program NS