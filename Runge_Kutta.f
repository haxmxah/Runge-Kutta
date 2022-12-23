c       MARTA XIULAN ARIBÓ HERRERA
C       PRÀCTICA 8 DE FÍSICA COMPUTACIONAL
C       En aquesta pràctica resoldrem l'equació d'Schrödinger independent 
c       del temps per trobar els autovalors i autovectors d'una partícula
c       mitjantçant el mètode de Runge Kutta 4
C-----------------------------------------------------------------------

        PROGRAM EQUACIO_SCHRODINGER
        Implicit none
c--------------------------- DEFINICIÓ DE VARIABLES --------------------
        integer i,c, pasos,neq, II
        parameter(pasos=500, neq=2)
        double precision E,beta,mew2,h2me,hw,xmin,xmax,L,dx,x,sigma
        parameter(mew2=1.d0)!eV/A^2
        parameter(h2me=3.80995d0)!eVA^2
        parameter(hw = 2.76041d0)!eV
        parameter (sigma = 1.175d0)
        common/energia/E
        common/parametre/beta
        double precision E1,E2,E3,E4,E5,E6,E_v(4),E11(3),E22(3),E33(3)
        double precision yyin(neq),phi(neq),phi1(neq),phi2(neq),
     *phi3(neq),Integral,autovalors(pasos),beta_v(3),dx2,Integral2
        open(1,file="auxiliar.dat")
        open(2,file="auxiliar2.dat")
        open(3,file="P8-20-21-b1-res1.dat")
c--------------------------- PRIMER APARTAT  ---------------------------
c       Obtenim les solucions corresponents per a l'equació diferencial amb
c       les condicions inicials establertes.

        Write(1,*)"#APARTAT 1"
        E1 = 0.43d0*hw !energies 
        E2 = 0.48d0*hw
        E3 = 1.40d0*hw
        E4 = 1.45d0*hw

        E_v = (/E1,E2,E3,E4/) !els posem dins d'un vector perquè sigui més compacte el codi

        L = 16.d0 !A
        xmin = -L/2.d0 !límits
        xmax = +L/2.d0

        dx = (xmax-xmin)/dble(pasos)!intervals

        beta= 0.d0!el passem per un common

        do c= 1, 4
            E = E_v(c)!el passem com un common
            Write(1,200)"Iteracio", "X", "Autovalor phi"
            Write(1,*)"#Solucions per l'energia E = ",c,E_v(c)
            yyin(1)= 0.d0 !condicions inicials
            yyin(2)= 10.d0**(-5.d0)

            do i = 1, pasos,1 !rRK amb 500 pasos
                x = xmin + dx*(i-1)
                call miRungeKutta4(x,dx,yyin,neq,phi)
                yyin = phi
                write(1,*) i,x,phi
            enddo
            write(1,"(a1)")
            write(1,"(a1)")
        enddo

        call system ("gnuplot -p plot1.gnu")!representació gràfica

c--------------------------- SEGON APARTAT  ---------------------------
c       Trobem els tres primers autovalors del sistema amb el mètode de 
c       la secant.
c-----------------------------------------------------------------------
        E5 = 2.d0*hw !energies
        E6 = 2.1d0*hw 

        E11 = (/E1,E3,E5/) !Utilitzem dos vectors de les energies per
        E22 = (/E2,E4,E6/) !fer més compacte el codi 

        DO c = 1, 3
            write(1,*)"#Autovalors per E =",c, E_v(c)
            write(2,*)"Convergència per l'autovalor",c !guardem a l'arxiu auxiliar 2
            DO II = 1,100 !Posem 100 per exemple

            yyin(1)= 0.d0 !condicions inicials
            yyin(2)= 10.d0**(-5.d0)

            do i = 1, pasos !obtenim phi1

                E = E11(C)!passem per un common l'energia
                x = xmin + dx*(i-1)
                call miRungeKutta4(x,dx,yyin,neq,phi1)
                yyin = phi1
     
            enddo

            yyin(1)= 0.d0!renovem les condicions inicials
            yyin(2)= 10.d0**(-5.d0)

            do i = 1, pasos !obtenim phi2

                E = E22(C)!passem per un common l'energia
                x = xmin + dx*(i-1)
                call miRungeKutta4(x,dx,yyin,neq,phi2)
                yyin = phi2
     
            enddo

            yyin(1)= 0.d0!renovem les condicions inicials
            yyin(2)= 10.d0**(-5.d0)

            E33(c)= (E11(c)*phi2(1)-E22(c)*phi1(1))/ (phi2(1)-phi1(1)) 
            !mètode de la secant
            Integral = 0.d0

            do i = 1, pasos !obtenim els autovectors i a l'hora calculem la seva
                            !normalització

                E = E33(c)!passem per un common l'energia
                x = xmin +dx*(i-1)
                call miRungeKutta4(x,dx,yyin,neq,phi3)!generem els autovalors phi3
                autovalors(i)=phi3(1) !guardem els valors dels autovalors
                yyin=phi3

                !integrem per normalitzar amb el mètode de trapezis
                if ((i.EQ.1).OR.(i.EQ.pasos)) then
                        Integral = Integral + (dx/2.d0)*(phi3(1)**2)
                else
                        Integral = Integral + (phi3(1)**2)*dx
                endif

            enddo
            write(2,*) II,E33(c),phi3(1) !guardem en un arxiu auxiliar per 
            !estudiar la convergència i gràficar

            if (abs(phi3(1)).LT.(10.d0**(-8.d0))) then !condició convergència
            write(1,"(a1,10x,a30)")"x","autovalors normalitzats"
            print*,"S'ha arribat a la convergència",II,phi3(1),
     *E33(C)
                do i = 1, pasos
                    x = xmin +dx*(i-1)
                    write(1,*)x,autovalors(i)/sqrt(integral) !escribim en l'arxiu 
                enddo
                write(1,"(a1)")
                write(1,"(a1)")
                exit

            endif

            E11(c) = E22(c) !sinó es cumpleixen les condicions
            E22(c) = E33(c) !tornem al bucle amb aquestes noves condicions

            ENDDO
            write(2,"(a1)")
            write(2,"(a1)")
        ENDDO

        call system ("gnuplot -p plot2.gnu")!representació gràfica
        call system ("gnuplot -p plot3.gnu")!representació gràfica

c--------------------------- TERCER APARTAT  ---------------------------
c       Càlculem de nou per l'estat fonamental els diferents autovalors 
c       i autovectors normalitzats per diferents betas
        
        L = 10.d0
        xmin = -L/2.d0 !límits
        xmax = +L/2.d0

        dx = (xmax-xmin)/dble(pasos)!intervals

        beta_v = (/0.d0,0.1d0,0.25d0/)!eV/A^4 
        !Els fiquem en un vector i els recorrem per compactar el codi

        E1 = 0.43d0*hw !prop de l'estat fonamental
        E2 = 0.48d0*hw

        write(3,*) "#Càlcul de probabilitat d'una particula" !titulillo
        DO c = 1,3
            beta = beta_v(c)

            DO II = 1, 100
            yyin(1)= 0.d0
            yyin(2)= 10.d0**(-5.d0)

            do i = 1, pasos !obtenim phi1

              E = E1!passem per un common l'energia
              x = xmin + dx*(i-1)
              call miRungeKutta4(x,dx,yyin,neq,phi1)
              yyin = phi1
     
            enddo

            yyin(1)= 0.d0
            yyin(2)= 10.d0**(-5.d0)

            do i = 1, pasos !obtenim phi1

               E = E2!passem per un common l'energia
               x = xmin + dx*(i-1)
               call miRungeKutta4(x,dx,yyin,neq,phi2)
               yyin = phi2
     
            enddo
            
            E3= (E1*phi2(1)-E2*phi1(1))/(phi2(1)-phi1(1)) !mètode de la secant

            yyin(1)= 0.d0
            yyin(2)= 10.d0**(-5.d0)
c            print*,E3
            Integral = 0.d0 
            do i = 1, pasos

                E = E3!passem per un common l'energia
                x = xmin +dx*(i-1)

                call miRungeKutta4(x,dx,yyin,neq,phi3)!generem els autovalors phi3
                autovalors(i)=phi3(1) !guardem els valors dels autovalors
                yyin=phi3

                !integrem per normalitzar amb el mètode de trapezis
                if ((i.EQ.1).OR.(i.EQ.pasos)) then
                        Integral = Integral + (dx/2.d0)*(phi3(1)**2)
                else
                        Integral = Integral + (phi3(1)**2)*dx
                endif

            enddo

            if (abs(phi3(1)).LT.(10.d0**(-8.d0))) then !condició convergència
            write(1,"(a1,10x,a30)")"x","autovalors normalitzats"
            print*,"S'ha arribat a la convergència",II,phi3(1),
     *E33(C)   
            dx2 = (2.d0*sigma)/pasos!aprofitem per calcular la probabilitat
            Integral2 = 0.d0 
c      Aprofitem el bucle per calcular la probabilitatt entre els intervals
c      [-sigma,sigma] amb les condicions de que x es trobi dins de l'interval
c      aquest, i amb el mètode de trapezis.

            do i = 1, pasos
                if ((x.GE.-sigma).AND.(x.LE.sigma))then
                    if ((x.GT.-sigma).AND.(i.LT.sigma)) then
                        Integral2= Integral2+ 
     *((autovalors(i)/sqrt(integral))**2)*dx2
                    else                        
                        Integral2= Integral2+ (dx2/2.d0)* 
     *((autovalors(i)/sqrt(integral))**2)
                    endif
                endif
                x = xmin +dx*(i-1)
                write(1,*)x,autovalors(i)/sqrt(integral)
            enddo
            print*, "La probabilitat per", beta, "és",integral2
            write(3,*)"La probabilitat per",beta,"és",integral2 !guardem les probabilitats
            write(1,"(a1)")
            write(1,"(a1)")
            exit
            endif
            E1 = E2
            E2 = E3

            ENDDO
        ENDDO
        call system ("gnuplot -p plot4.gnu")!gràfica dels autovectors segons beta




c--------------------------- FORMATS------------------------------------
200        format(a8,10x,a1,10x,a20)

        END PROGRAM 

C-----------------------------------------------------------------------

c--------------------------- SUBRUTINES I FUNCIONS ---------------------
C                Funcions utilitzades en aquesta prepràctica 
c-----------------------------------------------------------------------

           
c------------------- subrutina general de Runge Kutta ------------------
C       Té per INPUTS
c           |yyin = condicions inicials
c           |nequs = numero d'equacions
c           |t, dt = pas i temps
c           |yyout = resultat per un pas
c       Té per OUTPUTS 
c           |yyout = y1 (un pas de runge kutta)
C       SUBRUTINA TRETA DE BRUNO JULIÀ DEL CV 
c-----------------------------------------------------------------------

        subroutine miRungeKutta4(t,dt,yyin,nequs,yyout)
        implicit none
        integer nequs ,j
        double precision yyin(nequs), yyout(nequs),t,dt
        double precision k1(nequs),k2(nequs),k3(nequs),k4(nequs)
        double precision kk(nequs), EE
        common/energia/EE

        call derivada(t,yyin,K1,nequs)

        kk = yyin + (dt*k1/2.d0)
        call derivada(t+dt/2.d0,kk,k2,nequs)
    
        kk = yyin + (dt*k2/2.d0)
        call derivada(t+dt/2.d0,kk,k3,nequs)
        
        kk= yyin+ (dt*k3)
        call derivada(t+dt,kk,k4,nequs)

        yyout = yyin + dt*(1.d0/6.d0)*(k1+2.d0*k2+2.d0*k3+k4)

        end subroutine

c---------------------------SUBRUTINA DERIVADA--------------------------
c       !Aquesta subrutina és molt important!
c       Té per entrades yin on |yin(1) = angle 
c                              |yin(2) = derivada de l'angle 
c       Té per surtida yout on |dyout(1) = derivda de l'angle = yin(2)
c                              |dyout(2) = derivada de yin(2) = derivada
c                                                           segona de l'angle
c       En el nostre cas el nostre vector funció té dimensió 2
c       SUBRUTINA FETA A PARTIR DE LA IDEA DE BRUNO JULIÀ A LES CLASSES
C-----------------------------------------------------------------------

        subroutine derivada (t,yin,dyout,nequs) 
        implicit none
        integer nequs
        double precision t,yin(nequs),dyout(nequs),EE,bbeta,mew2,h2me,V
        common/energia/EE
        common/parametre/bbeta
        parameter(mew2=1.d0)!eV/A^2
        parameter(h2me=3.80995)!eVA^2

            dyout(1) = yin(2)
            V = (0.5*mew2)*(t**2.d0)+(bbeta*(t**4.d0))
            dyout(2) = (V-EE)*yin(1)*(1.d0/h2me)!terme

        end 






