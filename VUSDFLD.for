      subroutine vusdfld(nblock,nstatev,nfieldv,nprops,ndir,
     &	                 nshr,jElem,kIntPt,kLayer,kSecPt, 
     &                   stepTime,totalTime,dt,cmname, 
     &                   coordMp,direct,T,charLength,props, 
     &                    stateOld,stateNew,field )
C
C      include 'ABA_PARAM.INC'
C
      include 'vaba_param.inc'
C
      dimension jElem(nblock), coordMp(nblock,*), 
     1          direct(nblock,3,3), T(nblock,3,3), 
     2          charLength(nblock), props(nprops), 
     3          stateOld(nblock,nstatev), 
     4          stateNew(nblock,nstatev),
     5          field(nblock,nfieldv)
      character*80 cmname
      
      parameter( nrData=6 )
      character*3 cData(maxblk*nrData)
      dimension rData(maxblk*nrData), jData(maxblk*nrData)
C	  
	  real*8,dimension(6)::sigmaNew
	  real*8,dimension(1)::Critical_Strain	
	  real*8,dimension(1)::MaxStress	  
	  real*8,dimension(1)::term1
	  real*8,dimension(1)::term2
	  real*8,dimension(1)::term3
	  real*8,dimension(1)::term4
	  real*8,dimension(1)::term5
	  real*8,dimension(1)::maxprincipalE
	  real*8,dimension(1)::minprincipalE	
C	
	  jStatus = 1
	  call vgetvrm('LE',rdata,jdata,cdata,jStatus)
      if( jStatus .ne. 0 ) then
        call xplb_abqerr(-2,'Utility routine VGETVRM '//
     .      'failed to get variable.',0,zero,' ')
        call xplb_exit
      end if
C
	  Critical_Strain(1)=0.3000
C
      do k = 1, nblock
        StrainE11 = rData(k)
        StrainE22 = rData(nblock+k)	
        StrainE33 = rData(2*nblock+k)
        StrainE12 = rData(3*nblock+k)	

		term1(1) = (StrainE11+StrainE22)/2
		term2(1) = (StrainE11-StrainE22)/2
		term3(1) = StrainE12/2
		term4(1) = ((term2(1)*term2(1))+(term3(1)*term3(1)))
		term5(1) = SQRT(term4(1))
		maxprincipalE(1) = term1(1) + term5(1)
		minprincipalE(1) = term1(1) - term5(1)
		
        if (maxprincipalE(1) .ge. Critical_Strain(1)) then
         stateNew(k,1) = 0
        end if
      end do
      RETURN
      END