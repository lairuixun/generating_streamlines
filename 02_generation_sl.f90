!****************************************************
!	The program is used to generate streamline from velocity field.
!	Users can chenge their data with the domain of interest.
!	The code can be compiled using gfortran by:
!	gfortran sl.f90 -o sl 
!
!	lairuixun@163.com
!	Jan. 2020
!****************************************************

subroutine InterpolateVelocity(numNodes, nodesX,nodesY,nodesU,nodesV, searchX,searchY, radius,numPR,searchU,searchV)
	implicit none 
	integer			::i,j,k
	integer			::numNodes	
	real(kind=8),dimension(numNodes),intent(in)	::nodesX
	real(kind=8),dimension(numNodes),intent(in)	::nodesY
	real(kind=8),dimension(numNodes),intent(in)	::nodesU
	real(kind=8),dimension(numNodes),intent(in)	::nodesV
	real(kind=8),intent(in)						::searchX
	real(kind=8),intent(in)						::searchY
	real(kind=8),intent(out)					::searchU
	real(kind=8),intent(out)					::searchV
	real,intent(in) 							::radius 
	real			 							::distance
	integer,intent(out)							::numPR	! number of points within the circle 
	
	!searching nodes within the circle 
	!and caculate the velocity of the searching point 
	do i=1, numNodes
		distance = ((searchX-nodesX(i))**2.0 + (searchY-nodesY(i))**2.0 )**0.5
		! caculate the velocity of the searching points
		if (distance .lt. radius) then 
			numPR = numPR+1
			searchU	= (searchU + nodesU(i))/numPR 
			searchV	= (searchV + nodesV(i))/numPR 
			
			!print *, nodesX(i),nodesY(i),searchU,searchV,numPR 
		
		end if 
		
	end do 
	if (numPR == 0) then 
		print *, "there's no nodes searched in the circle; Caculating the velocity of the searching point failed."
	end if 

end subroutine InterpolateVelocity

program generate_sl

	implicit none

	character(50)	::tempstring
	integer			::i,j,k
	integer			::numNodes
	real(kind=8),allocatable,dimension(:)::nodesX
	real(kind=8),allocatable,dimension(:)::nodesY
	real(kind=8),allocatable,dimension(:)::nodesU
	real(kind=8),allocatable,dimension(:)::nodesV
	real(kind=8)			::searchX	!x coordinate of seed point 
	real(kind=8)			::searchY	!y coordinate of seed point 
	real(kind=8)			::searchU
	real(kind=8)			::searchV
	real 					::radius 
	real 					::distance
	integer					::numPR	! number of points within the circle 
	real 					::stepsize
	real					::DeltaR
	!boundary of the domain
	real 					::Domain_minX,Domain_maxX,Domain_minY,Domain_maxY
	real(kind=8)			::trajX,trajX2
	real(kind=8)			::trajY,trajY2
	real(kind=8)			::N1,N2,N3,N4
	real(kind=8)			::M1,M2,M3,M4
	
	!allocate dimension of nodes
	! 7216 is the number of nodes 
	numNodes 		= 7216
	allocate(nodesX(numNodes))
	allocate(nodesY(numNodes))
	allocate(nodesU(numNodes))
	allocate(nodesV(numNodes))
	Domain_minX =473262
	Domain_maxX =477555
	Domain_minY =4174273
	Domain_maxY =4178678
	searchX		=477249.0
	searchY		=4178183.0
	searchU 	= 0.0
	searchV 	= 0.0
	radius 		= 50.0
	numPR 		= 0
	stepsize	= 0.2
	DeltaR		= 100.0
	trajX = searchX
	trajY = searchY
	
	!read nodes and velocity values
	open(1,file='01_velocity_field.plt',status='old')
		read(1,*) tempstring
		read(1,*) tempstring
		read(1,*) tempstring
		do i=1,numNodes
			read(1,*) nodesX(i),nodesY(i),nodesU(i),nodesV(i)
			!print *, nodesX(i),nodesY(i)
		end do
	close(1)

	open(2,file='traj.txt',status='old')
	
	!do while acquiring enough streamline point 
	!interpolate velocity using 4th order Runge-Kutta
	j = 0
	do   while (j .le. 20000)
		!N = DeltaR * v /sqrt(u**2+v**2)
		!M = DeltaR * u /sqrt(u**2+v**2)
		call InterpolateVelocity(numNodes, nodesX,nodesY,nodesU,nodesV, trajX,trajY, radius,numPR,searchU,searchV)	
		N1 = DeltaR*searchV/(searchU**2.0+searchV**2.0)**0.5
		M1 = DeltaR*searchU/(searchU**2.0+searchV**2.0)**0.5
		searchX = searchX + 0.5*stepsize*M1
		searchY = searchY + 0.5*stepsize*N1
		call InterpolateVelocity(numNodes, nodesX,nodesY,nodesU,nodesV, trajX,trajY, radius,numPR,searchU,searchV)		
		N2 = (DeltaR+0.5*stepsize)*searchV/(searchU**2.0+searchV**2.0)**0.5
		M2 = (DeltaR+0.5*stepsize)*searchU/(searchU**2.0+searchV**2.0)**0.5
		searchX = searchX + 0.5*stepsize*M2
		searchY = searchY + 0.5*stepsize*N2
		call InterpolateVelocity(numNodes, nodesX,nodesY,nodesU,nodesV, trajX,trajY, radius,numPR,searchU,searchV)		
		N3 = (DeltaR+0.5*stepsize)*searchV/(searchU**2.0+searchV**2.0)**0.5
		M3 = (DeltaR+0.5*stepsize)*searchU/(searchU**2.0+searchV**2.0)**0.5
		searchX = searchX + 0.5*stepsize*M3
		searchY = searchY + 0.5*stepsize*N3
		call InterpolateVelocity(numNodes, nodesX,nodesY,nodesU,nodesV, trajX,trajY, radius,numPR,searchU,searchV)		
		N4 = (DeltaR+0.5*stepsize)*searchV/(searchU**2.0+searchV**2.0)**0.5
		M4 = (DeltaR+0.5*stepsize)*searchU/(searchU**2.0+searchV**2.0)**0.5
		trajY2 = trajY + stepsize*(N1+2.0*N2+2.0*N3+N4)/6.0
		trajX2 = trajX + stepsize*(M1+2.0*M2+2.0*M3+M4)/6.0
		
		write(2,*) trajX2,trajY2 
		j = j + 1		
		print *, j, trajX2,trajY2 

		trajX = trajX2
		trajY = trajY2
		
		if ((trajX .lt. Domain_minX) .or.  (trajX .gt. Domain_maxX) ) then 
			exit
		end if 
		if ((trajY .lt. Domain_minY) .or.  (trajY .gt. Domain_maxY) ) then 
			exit
		end if 		
	end do
	
	close(2)
	
	
end program
