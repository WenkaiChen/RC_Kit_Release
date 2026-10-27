      subroutine align_rimers (coord1,coord2,n_atom,n_mode,atmnu1,atmnu2
     &,l_matrix, logid)



      implicit real*8 (a-h,o-z)
      real*8 coord1(n_mode),coord2(n_mode)
      integer*4 atmnu1(n_atom),atmnu2(n_atom)
      real*8 x1(3,n_atom),x2(3,n_atom),c2(3,n_atom)
      real*8 maxerr,maxemin
      real*8, allocatable :: x2a(:,:),wk(:,:), temp_l(:,:)
      real*8 l_matrix(n_mode,n_mode)
      integer*4 ord(n_atom),swap(3,6)
      integer*4 i,j, jj, iforder, error, logid
      logical loconly
      data swap / 1,2,3, 1,3,2, 2,1,3, 2,3,1, 3,1,2, 3,2,1/
      allocate(x2a(3,n_atom))
      allocate(wk(n_atom,n_atom))
      allocate(temp_l(n_mode,n_mode))
      iforder = 2

      open(2333,file="ALIGN_GEOM.TXT",status="REPLACE")
      write (2333,*) 'reordering atoms'

      J=0
      DO I=1,N_ATOM,1
        X1(1,I)=COORD1(J+1)
        X1(2,I)=COORD1(J+2)
        X1(3,I)=COORD1(J+3)
        x2(1,I)=COORD2(J+1)
        x2(2,I)=COORD2(J+2)
        x2(3,I)=COORD2(J+3)
      J=J+3
      END DO
      do i= 1,n_atom
	x2a(1,i)= x2(1,i)
	x2a(2,i)= x2(2,i)
	x2a(3,i)= x2(3,i)

      end do

      if (iforder.eq.0) then

	do i= 1,n_atom
	  ord(i)= i
          write (2333,*) iforder
	  write (2333,*) i,ord(i)
	end do



      end if

      if (iforder.ne.2) then
        loconly= .true.
        maxemin= 1.D30
      else
        loconly= .false.
        maxemin= 0.D0
      end if
      nswap= 0



      write (2333,*) 'initial cartes coords mol 1'
      do i= 1,n_atom
	write (2333,'(i4,3f10.5)') atmnu1(i),(x1(j,i),j=1,3)
      end do
      write (2333,*) 'initial cartes coords mol 2'
      do i= 1,n_atom
	write (2333,'(i4,3f10.5)') atmnu2(i),(x2(j,i),j=1,3)
      end do





77    continue
      nswap= nswap + 1


      kmap= -1
78    continue
      kmap= kmap + 1

79    continue
      write (2333,'(a,2i2)') ' loop',nswap,kmap

      iswx= 1
      if (mod(kmap,2).eq.1) iswx= -1
      iswy= 1
      if (mod(kmap,4).eq.2.or.mod(kmap,4).eq.3) iswy= -1
      iswz= 1
      if (mod(kmap,8)/4.eq.1) iswz= -1
      iswap= swap(1,nswap)
      jswap= swap(2,nswap)
      kswap= swap(3,nswap)
      maxerr= 0.d0

      do k= 1,n_atom
	x2(iswap,k)= x2a(1,k) * iswx
	x2(jswap,k)= x2a(2,k) * iswy
	x2(kswap,k)= x2a(3,k) * iswz
      end do

      write (2333,*) 'mod cartes coords mol 2'
      do i= 1,n_atom
	write (2333,'(i4,3f10.5)') atmnu2(i),(x2(k,i),k=1,3)
      end do

      if (loconly) then


	maxerr= 0.d0
        do i= 1,n_atom
	  j= ord(i)
	  xxx= (x2(1,i)-x1(1,j))**2 + (x2(2,i)-x1(2,j))**2 +
     $	       (x2(3,i)-x1(3,j))**2
	  if (maxerr.lt.xxx) then
	    maxerr= xxx
	    jmax= j
	  end if
	end do
	maxerr= sqrt(maxerr)
	write (2333,'(a,f6.2,a,i4)') ' read in order, maxerr=',maxerr,
     $		' for atom',jmax
	if (maxerr.lt.maxemin) then
	  maxemin= maxerr
	  kmapmin= kmap
	  nswapmin= nswap
	end if

      else



        do i= 1,n_atom
	  ord(i)= 0
        end do
	if (nswap.eq.1 .and. kmap.eq.0) write (2333,*) iswx,iswy,iswz

        do i= 1,n_atom
	  do j= 1,n_atom
	    xxx= (x2(1,i)-x1(1,j))**2 + (x2(2,i)-x1(2,j))**2 +
     $	         (x2(3,i)-x1(3,j))**2
	    if (atmnu1(j).eq.atmnu2(i) .and. xxx.lt.2) ord(i)= j
          end do
        end do

      end if

      do i= 1,n_atom
	if (maxemin.gt.0.3 .or. ord(i).eq.0) then
	  if (ord(i).eq.0) write (2333,*) 'UNMAPPED #2 atom ',i,atmnu2(i)


	  if (kmap.lt.7) goto 78


	  if (nswap.lt.6) goto 77



	    if (loconly) then
	      if (maxemin.gt.30.) then
	        do j= 1,n_atom
	          write (2333,'(i3,3x,3f10.6,3x,i3,3f10.6)') atmnu1(j),
     $		  (x1(k,j),k=1,3), atmnu2(j),(x2(k,j),k=1,3)
	        end do
	        stop 'Unmapped atom'
	      else

		kmap= kmapmin
		nswap= nswapmin
		write (2333,*) 'settling for map=',kmap,' swap=',kswap,
     $			' maxerr=',maxemin
            if (maxemin.gt.30.) then
              write (2333,*) 'Warning: Geom Displacment VERY BIG!!!'
              write (LOGID,*) 'settling for map=',kmap,' swap=',kswap,
     $			' maxerr=',maxemin
              write (LOGID,*) 'Warning: Geom Displacment VERY BIG!!!'
            end if
		maxemin= 0.0
		goto 79
	      end if

	    else

              WRITE(LOGID,*) 'The geometry displacement is large,'
              WRITE(LOGID,*) 'the program cannot rotate to a $
     $reasonable configuration.'
              WRITE(LOGID,*) 'The atom order will be read from $
     $a file named "align_order".'
              WRITE(2333,*) 'The geometry displacement is quite $
     $large, '
              WRITE(2333,*) 'the program cannot rotate to a $
     $reasonable configuration.'
              WRITE(2333,*) 'The atom order will be read from $
     $a file named "align_order".'
 	      write (2333,*) 'TRYING LOCATION ONLY--- READ IN ORDER'
              open(5,file="align_order",STATUS="OLD",IOSTAT=ERROR)
              IF (ERROR/=0) THEN

                STOP "File 'align_order' does not exist, exit."
              END IF
	      read (5,*) (ord(j),j=1,n_atom)
              close(5)
	      write (2333,*) 'input order=',(ord(j),j=1,n_atom)
	      iord= 0
	      do j= 1,n_atom
		j1= j
		if (atmnu1(j1).ne.atmnu2(j)) iord= 0
	      end do
	        if (iord.eq.0) then
		  write (2333,*) 'ORDER ERROR in atomic numbers'
	          do jj= 1,n_atom
	            write (2333,'(i3,3x,3f10.6,3x,i3,3f10.6)') atmnu1(jj),
     $		    (x1(k,jj),k=1,3), atmnu2(jj),(x2(k,jj),k=1,3)
	          end do

		end if
	      loconly= .true.
	      maxemin= 1.e30
	      nswap= 0
	      goto 77
	    end if

	end if
      end do

      write (2333,*) 'final axis swaps:',iswx,iswy,iswz
      write (2333,*) 'final xyz swap pattern=',nswap
      write (2333,*) 'reorder list:'
      write (2333,'(20i4)') (ord(i),i=1,n_atom)



      do i= 1,n_atom
	i1= ord(i)
	wk(1,i1)= x2(1,i)
	wk(2,i1)= x2(2,i)
	wk(3,i1)= x2(3,i)
      end do

      do i= 1,n_atom
	do j= 1,3
	  x2(j,i)= wk(j,i)
	end do
      end do

      do i= 1,n_atom
	i1= ord(i)
	wk(1,i1)= atmnu2(i)
      end do

      do i= 1,n_atom
	atmnu2(i)= nint(wk(1,i))
      end do



      if (n_mode.eq.0) return

      write (2333,*) 'axis signs=',iswx,iswy,iswz


      JJ=0
      DO J=1,N_MODE,1
        DO I=1,N_ATOM,1
          JJ=ORD(I)
          TEMP_L(3*JJ-3+ISWAP,J) = L_MATRIX(3*I-2,J)*iswx
          TEMP_L(3*JJ-3+JSWAP,J) = L_MATRIX(3*I-1,J)*iswy
          TEMP_L(3*JJ-3+KSWAP,J) = L_MATRIX(3*I  ,J)*iswz
        END DO
      END DO
      L_MATRIX=TEMP_L




      WRITE(2333,"('COORDS1 AFTER TRANSFORMATION:')")
      DO I=1,N_ATOM,1
        WRITE(2333,"(I4,3F10.4)") ATMNU1(I), (x1(J,I),J=1,3)
      END DO
      WRITE(2333,"('COORDS2 AFTER TRANSFORMATION:')")
      DO I=1,N_ATOM,1
        WRITE(2333,"(I4,3F10.4)") ATMNU2(I), (x2(J,I),J=1,3)
      END DO
      J=0
      DO I=1,N_ATOM,1
        COORD2(J+1)=x2(1,I)
        COORD2(J+2)=x2(2,I)
        COORD2(J+3)=x2(3,I)
        J=J+3
      END DO
      WRITE(2333,"('TEST COORD2:')")
      WRITE(2333,"(3F10.4)") COORD2
      deallocate (x2a)
      deallocate (wk)
      deallocate (temp_l)
      close(2333)
      return
      end
