!**********************************************************************************************
! Subroutines for JACOBI iteration for the Crank--Nicholson and pressure steps.
!
! These subroutines use flux-conservative differencing, and the cell-averaged viscosities are 
! copied from the master subroutine.
!
!**********************************************************************************************

subroutine do_jacobi_u(u3,u3_old,RHS,viscosity,density,dx,dy,dz,dt,ex,ey,sx,sy,maxl,maxn,ex_max)
  implicit none
  
  integer :: sx,sy,ex,ey,flg
  integer :: i,j,k,maxl,maxn,iteration_time,ex_max,ex_loc
  
  double precision :: u3(sx-1:ex+1,sy-1:ey+1,0:maxn-2), u3_old(sx-1:ex+1,sy-1:ey+1,0:maxn-2), &
       RHS(sx-1:ex+1,sy-1:ey+1,0:maxn-2)
  double precision :: viscosity(sx-2:ex+2,sy-2:ey+2,0:maxn)
  double precision :: density(sx-1:ex+1,sy-1:ey+1,0:maxn)
  double precision :: dx,dy,dz,dt,diag_val,residual,ax,ay,az
  double precision :: mu_plushalf_x_val,mu_minushalf_x_val,mu_plushalf_y_val,mu_minushalf_y_val, &
       mu_plushalf_z_val,mu_minushalf_z_val
  
  double precision :: u_minusz,u_plusz,rho_ugrid
  
  ax=dt/(dx*dx)
  ay=dt/(dy*dy)
  az=dt/(dz*dz)
  
  if(ex==ex_max)then
     ex_loc=ex-1
  else
     ex_loc=ex
  end if
  
  do k=0,maxn-2
     do j=sy,ey
        do i=sx,ex_loc
           
           if(k.eq.0)then
              u_minusz=-2.d0*u3_old(i,j,0)+u3_old(i,j,1)/3.d0
           else
              u_minusz=u3_old(i,j,k-1)
           end if
           
           if(k.eq.(maxn-2))then
              u_plusz=-2.d0*u3_old(i,j,maxn-2)+u3_old(i,j,maxn-3)/3.d0
           else
              u_plusz=u3_old(i,j,k+1)
           end if
           
           mu_plushalf_x_val =viscosity(i+1,j+1,k+1)
           mu_minushalf_x_val=viscosity(i,  j+1,k+1)
           
           mu_plushalf_y_val =(viscosity(i,j+1,k+1)+viscosity(i+1,j+1,k+1)+viscosity(i,j+2,k+1)+viscosity(i+1,j+2,k+1))/4.d0
           mu_minushalf_y_val=(viscosity(i,j+1,k+1)+viscosity(i+1,j+1,k+1)+viscosity(i,j  ,k+1)+viscosity(i+1,j  ,k+1))/4.d0
           
           mu_plushalf_z_val =(viscosity(i,j+1,k+1)+viscosity(i+1,j+1,k+1)+viscosity(i,j+1,k+2)+viscosity(i+1,j+1,k+2))/4.d0
           mu_minushalf_z_val=(viscosity(i,j+1,k+1)+viscosity(i+1,j+1,k+1)+viscosity(i,j+1,k  )+viscosity(i+1,j+1,k  ))/4.d0
           
           rho_ugrid=(density(i+1,j+1,k+1)+density(i,j+1,k+1))/2.d0
           
           diag_val=1.d0+(ax/2.d0)*(1.d0/rho_ugrid)*(mu_plushalf_x_val+mu_minushalf_x_val)+&
                (ay/2.d0)*(1.d0/rho_ugrid)*(mu_plushalf_y_val+mu_minushalf_y_val)+&
                (az/2.d0)*(1.d0/rho_ugrid)*(mu_plushalf_z_val+mu_minushalf_z_val)
           
           residual = (1.d0/diag_val) * ( &
                (ax/2.d0)*(1.d0/rho_ugrid)*(mu_plushalf_x_val*u3_old(i+1,j,k)+mu_minushalf_x_val*u3_old(i-1,j,k)) + &
                (ay/2.d0)*(1.d0/rho_ugrid)*(mu_plushalf_y_val*u3_old(i,j+1,k)+mu_minushalf_y_val*u3_old(i,j-1,k)) + &
                (az/2.d0)*(1.d0/rho_ugrid)*(mu_plushalf_z_val*u_plusz    +mu_minushalf_z_val*u_minusz   ) + RHS(i,j,k))
           
           
           u3(i,j,k) = residual
           
        end do
     end do
  end do
  
  return
end subroutine do_jacobi_u

! **********************************************************************************************

subroutine do_jacobi_v(v3,v3_old,RHS,viscosity,density,dx,dy,dz,dt,ex,ey,sx,sy,maxl,maxn,ex_max)
  implicit none
  
  integer :: sx,sy,ex,ey
  integer :: i,j,k,maxl,maxn,ex_max
  
  double precision :: v3(sx-1:ex+1,sy-1:ey+1,0:maxn-2), v3_old(sx-1:ex+1,sy-1:ey+1,0:maxn-2), &
       RHS(sx-1:ex+1,sy-1:ey+1,0:maxn-2),viscosity(sx-2:ex+2,sy-2:ey+2,0:maxn),             &
       density(sx-1:ex+1,sy-1:ey+1,0:maxn)
  double precision :: dx,dy,dz,dt,diag_val,residual,ax,ay,az
  double precision :: mu_plushalf_x_val,mu_minushalf_x_val,mu_plushalf_y_val,mu_minushalf_y_val, &
       mu_plushalf_z_val,mu_minushalf_z_val
  
  double precision :: v_minusz,v_plusz,v_plusx,v_minusx,rho_vgrid
  
  ax=dt/(dx*dx)
  ay=dt/(dy*dy)
  az=dt/(dz*dz)
  
  do k=0,maxn-2
     do j=sy,ey
        do i=sx,ex
           
           if(k.eq.0)then
              v_minusz=-2.d0*v3_old(i,j,0)+v3_old(i,j,1)/3.d0
           else
              v_minusz=v3_old(i,j,k-1)
           end if
           
           if(k.eq.(maxn-2))then
              v_plusz=-2.d0*v3_old(i,j,maxn-2)+v3_old(i,j,maxn-3)/3.d0
           else
              v_plusz=v3_old(i,j,k+1)
           end if
           
           v_plusx=v3_old(i+1,j,k)
           v_minusx=v3_old(i-1,j,k)
           
           mu_plushalf_x_val= (viscosity(i,j+1,k+1)+viscosity(i,j,k+1)+viscosity(i+1,j+1,k+1)+viscosity(i+1,j,k+1))/4.d0
           mu_minushalf_x_val=(viscosity(i,j+1,k+1)+viscosity(i,j,k+1)+viscosity(i-1,j+1,k+1)+viscosity(i-1,j,k+1))/4.d0
           
           mu_plushalf_y_val= viscosity(i,j+1,k+1)
           mu_minushalf_y_val=viscosity(i,j,  k+1)
           
           mu_plushalf_z_val= (viscosity(i,j+1,k+1)+viscosity(i,j,k+1)+viscosity(i,j+1,k+2)+viscosity(i,j,k+2))/4.d0
           mu_minushalf_z_val=(viscosity(i,j+1,k+1)+viscosity(i,j,k+1)+viscosity(i,j+1,k  )+viscosity(i,j,k  ))/4.d0
           
           rho_vgrid=(density(i,j+1,k+1)+density(i,j,  k+1))/2.d0
           
           diag_val=1.d0+(ax/2.d0)*(1.d0/rho_vgrid)*(mu_plushalf_x_val+mu_minushalf_x_val)+&
                (ay/2.d0)*(1.d0/rho_vgrid)*(mu_plushalf_y_val+mu_minushalf_y_val)+&
                (az/2.d0)*(1.d0/rho_vgrid)*(mu_plushalf_z_val+mu_minushalf_z_val)
           
           residual = (1.d0/diag_val) * ( &
                (ax/2.d0)*(1.d0/rho_vgrid)*(mu_plushalf_x_val*v3_old(i+1,j,k)+mu_minushalf_x_val*v3_old(i-1,j,k)) + &
                (ay/2.d0)*(1.d0/rho_vgrid)*(mu_plushalf_y_val*v3_old(i,j+1,k)+mu_minushalf_y_val*v3_old(i,j-1,k)) + &
                (az/2.d0)*(1.d0/rho_vgrid)*(mu_plushalf_z_val*v_plusz    +mu_minushalf_z_val*v_minusz   ) + RHS(i,j,k))
           
           v3(i,j,k) = residual
           
        end do
     end do
  end do
  
  return
end subroutine do_jacobi_v

! **********************************************************************************************

subroutine do_jacobi_w(w3,w3_old,RHS,viscosity,density,dx,dy,dz,dt,ex,ey,sx,sy,maxl,maxn,ex_max)
  implicit none
  
  integer :: sx,sy,ex,ey,flg
  integer :: i,j,k,maxl,maxn,iteration_time
  integer :: ex_max
  
  double precision :: w3(sx-1:ex+1,sy-1:ey+1,0:maxn-1), w3_old(sx-1:ex+1,sy-1:ey+1,0:maxn-1), &
       RHS(sx-1:ex+1,sy-1:ey+1,0:maxn-1),viscosity(sx-2:ex+2,sy-2:ey+2,0:maxn), &
       density(sx-1:ex+1,sy-1:ey+1,0:maxn)
  double precision :: dx,dy,dz,dt,diag_val,residual,ax,ay,az, &
       mu_plushalf_x_val,mu_minushalf_x_val,mu_plushalf_y_val,mu_minushalf_y_val, &
       mu_plushalf_z_val,mu_minushalf_z_val
  double precision :: w_plusx,w_minusx,rho_wgrid
  
  ax=dt/(dx*dx)
  ay=dt/(dy*dy)
  az=dt/(dz*dz)
  
  do k=1,maxn-2
     do j=sy,ey
        do i=sx,ex
           
           mu_plushalf_x_val= (viscosity(i,j+1,k)+viscosity(i,j+1,k+1)+viscosity(i+1,j+1,k)+viscosity(i+1,j+1,k+1))/4.d0
           mu_minushalf_x_val=(viscosity(i,j+1,k)+viscosity(i,j+1,k+1)+viscosity(i-1,j+1,k)+viscosity(i-1,j+1,k+1))/4.d0
           
           mu_plushalf_y_val= (viscosity(i,j+1,k)+viscosity(i,j+1,k+1)+viscosity(i,j+2,k)+viscosity(i,j+2,k+1))/4.d0
           mu_minushalf_y_val=(viscosity(i,j+1,k)+viscosity(i,j+1,k+1)+viscosity(i,j  ,k)+viscosity(i,j  ,k+1))/4.d0
           
           mu_plushalf_z_val=  viscosity(i,j+1,k+1)
           mu_minushalf_z_val= viscosity(i,j+1,k  )
           
           w_plusx=w3_old(i+1,j,k)
           w_minusx=w3_old(i-1,j,k)
           
           rho_wgrid=(density(i,j+1,k+1)+density(i,j+1,k))/2.d0
           
           diag_val=1.d0+(ax/2.d0)*(1.d0/rho_wgrid)*(mu_plushalf_x_val+mu_minushalf_x_val)+&
                (ay/2.d0)*(1.d0/rho_wgrid)*(mu_plushalf_y_val+mu_minushalf_y_val)+&
                (az/2.d0)*(1.d0/rho_wgrid)*(mu_plushalf_z_val+mu_minushalf_z_val)
           
           residual=(1.d0/diag_val)*( &
                (ax/2.d0)*(1.d0/rho_wgrid)*(mu_plushalf_x_val*w3_old(i+1,j,k)+mu_minushalf_x_val*w3_old(i-1,j,k)) + &
                (ay/2.d0)*(1.d0/rho_wgrid)*(mu_plushalf_y_val*w3_old(i,j+1,k)+mu_minushalf_y_val*w3_old(i,j-1,k)) + &
                (az/2.d0)*(1.d0/rho_wgrid)*(mu_plushalf_z_val*w3_old(i,j,k+1)+mu_minushalf_z_val*w3_old(i,j,k-1)) + RHS(i,j,k))
           
           w3(i,j,k) = residual
        end do
     end do
  end do
  
  return
end subroutine do_jacobi_w

! **********************************************************************************************

subroutine do_jacobi_p(pres,pres_temp,RHS,density,dx,dy,dz,dt,ex,ey,sx,sy,maxl,maxn,ex_max)
  implicit none
  
  integer :: sx,sy,ex,ey,flg
  integer :: i,j,k,maxl,maxn,iteration_time
  integer :: ex_max,sx_loc
  
  double precision :: pres(sx-1:ex+1,sy-1:ey+1,0:maxn), pres_temp(sx-1:ex+1,sy-1:ey+1,0:maxn), &
       RHS(sx-1:ex+1,sy-1:ey+1,0:maxn), density(sx-1:ex+1,sy-1:ey+1,0:maxn)
  double precision :: dx,dy,dz,dt,diag_val,residual,ax,ay,az
  double precision :: rho_plusx,rho_minusx,rho_plusy,rho_minusy,rho_plusz,rho_minusz
  
  ax=1.d0/(dx*dx)
  ay=1.d0/(dy*dy)
  az=1.d0/(dz*dz)
  
  do k=1,maxn-1
     do j=sy,ey
        do i=sx,ex
           
           rho_plusx= (density(i+1,j,k)+density(i,j,k))/2.d0
           rho_minusx=(density(i-1,j,k)+density(i,j,k))/2.d0
           
           rho_plusy= (density(i,j+1,k)+density(i,j,k))/2.d0
           rho_minusy=(density(i,j-1,k)+density(i,j,k))/2.d0
           
           rho_plusz= (density(i,j,k+1)+density(i,j,k))/2.d0
           rho_minusz=(density(i,j,k-1)+density(i,j,k))/2.d0
           
           diag_val=ax*( (1.d0/rho_plusx)+(1.d0/rho_minusx)) + &
                ay*( (1.d0/rho_plusy)+(1.d0/rho_minusy)) + &
                az*( (1.d0/rho_plusz)+(1.d0/rho_minusz))
           
           residual = (1.d0/diag_val) * ( &
                ay*( (1.d0/rho_minusy)*pres_temp(i,j-1,k) + (1.d0/rho_plusy)*pres_temp(i,j+1,k))  + &
                ax*( (1.d0/rho_minusx)*pres_temp(i-1,j,k) + (1.d0/rho_plusx)*pres_temp(i+1,j,k))  + &
                az*( (1.d0/rho_minusz)*pres_temp(i,j,k-1) + (1.d0/rho_plusz)*pres_temp(i,j,k+1))  - RHS(i,j,k))
           
           pres(i,j,k) = residual
        end do
     end do
  end do
  
  return
end subroutine do_jacobi_p

! **********************************************************************************************

subroutine get_difference(ua,ub,ex,ey,ez,sx,sy,sz,diff)
  implicit none
  
  double precision :: diff,sum
  integer :: ex,ey,ez,sx,sy,sz,i,j,k
  double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz:ez) :: ua, ub
  
  sum = 0.0D0
  
  do k=sz,ez
     do j = sy,ey
        do i = sx,ex
           sum = sum + (ua(i,j,k)-ub(i,j,k))**2
        end do
     end do
  end do
  
  diff=sum  
  
  Return
End subroutine get_difference

! **********************************************************************************************
