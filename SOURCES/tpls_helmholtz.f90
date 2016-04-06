module tpls_helmholtz

  implicit none

contains

  !**********************************************************************************************
  ! Subroutines for SOR iteration for the Crank--Nicholson step.
  !
  ! These subroutines use flux-conservative differencing, and the cell-averaged viscosities are 
  ! copied from the master subroutine.
  !
  !**********************************************************************************************
  
  subroutine do_sor_u(u3,RHS,viscosity,density,flg,iteration_time)

    use tpls_constants
    use tpls_mpi
    use mpi
    implicit none

    !----- Input -----!

    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn-2), intent(inout) :: u3
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn-2), intent(in)    :: RHS
    double precision, dimension(sx-2:ex+2,sy-2:ey+2,0:maxn)  , intent(in)    :: viscosity
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,0:maxn)  , intent(in)    :: density

    integer :: flg, iteration_time

    !----- Miscellaneous -----!

    integer :: i, j, k, ex_loc
    double precision :: diag_val, relax, residual, ax, ay, az
    double precision :: mu_plushalf_x_val, mu_minushalf_x_val
    double precision :: mu_plushalf_y_val, mu_minushalf_y_val
    double precision :: mu_plushalf_z_val,mu_minushalf_z_val
    double precision :: u_minusz,u_plusz,rho_ugrid
    
    !----- SOR iteration -----!
    
    relax = 1.2D+00
    
    ax = dt/(dx*dx)
    ay = dt/(dy*dy)
    az = dt/(dz*dz)
    
    if(ex == ex_max)then
       ex_loc = ex-1
    else
       ex_loc = ex
    end if
    
    if(mod(iteration_time,2).eq.0)then   
       
       do k = 0,maxn-2
          do j = sy,ey
             do i = mod(k+j+flg,2)+sx,ex_loc,2
                
                if(k.eq.0)then
                   u_minusz = -2.d0*u3(i,j,0)+u3(i,j,1)/3.d0
                else
                   u_minusz = u3(i,j,k-1)
                end if
                
                if(k.eq.(maxn-2))then
                   u_plusz = -2.d0*u3(i,j,maxn-2)+u3(i,j,maxn-3)/3.d0
                else
                   u_plusz = u3(i,j,k+1)
                end if
                
                mu_plushalf_x_val  = viscosity(i+1,j+1,k+1)
                mu_minushalf_x_val = viscosity(i,  j+1,k+1)
                
                mu_plushalf_y_val  = (viscosity(i,j+1,k+1)+viscosity(i+1,j+1,k+1)+viscosity(i,j+2,k+1)+viscosity(i+1,j+2,k+1))/4.d0
                mu_minushalf_y_val = (viscosity(i,j+1,k+1)+viscosity(i+1,j+1,k+1)+viscosity(i,j  ,k+1)+viscosity(i+1,j  ,k+1))/4.d0
                
                mu_plushalf_z_val  = (viscosity(i,j+1,k+1)+viscosity(i+1,j+1,k+1)+viscosity(i,j+1,k+2)+viscosity(i+1,j+1,k+2))/4.d0
                mu_minushalf_z_val = (viscosity(i,j+1,k+1)+viscosity(i+1,j+1,k+1)+viscosity(i,j+1,k  )+viscosity(i+1,j+1,k  ))/4.d0
                
                rho_ugrid = (density(i+1,j+1,k+1)+density(i,j+1,k+1))/2.d0
                
                diag_val = 1.d0+(ax/2.d0)*(1.d0/rho_ugrid)*(mu_plushalf_x_val+mu_minushalf_x_val) + &
                     (ay/2.d0)*(1.d0/rho_ugrid)*(mu_plushalf_y_val+mu_minushalf_y_val) + &
                     (az/2.d0)*(1.d0/rho_ugrid)*(mu_plushalf_z_val+mu_minushalf_z_val)
                
                residual = (1.d0/diag_val) * ( &
                     (ax/2.d0)*(1.d0/rho_ugrid)*(mu_plushalf_x_val*u3(i+1,j,k)+mu_minushalf_x_val*u3(i-1,j,k)) + &
                     (ay/2.d0)*(1.d0/rho_ugrid)*(mu_plushalf_y_val*u3(i,j+1,k)+mu_minushalf_y_val*u3(i,j-1,k)) + &
                     (az/2.d0)*(1.d0/rho_ugrid)*(mu_plushalf_z_val*u_plusz    +mu_minushalf_z_val*u_minusz   ) + RHS(i,j,k))
                
                u3(i,j,k) = u3(i,j,k)+relax*(residual-u3(i,j,k))
                
             end do
          end do
       end do
       
    else
       
       do k = 0,maxn-2
          do j = sy,ey
             do i = mod(j+flg,2)+sx,ex_loc,2
                
                if(k.eq.0)then
                   u_minusz = -2.d0*u3(i,j,0)+u3(i,j,1)/3.d0
                else
                   u_minusz = u3(i,j,k-1)
                end if
                
                if(k.eq.(maxn-2))then
                   u_plusz = -2.d0*u3(i,j,maxn-2)+u3(i,j,maxn-3)/3.d0
                else
                   u_plusz = u3(i,j,k+1)
                end if
                
                mu_plushalf_x_val  = viscosity(i+1,j+1,k+1)
                mu_minushalf_x_val = viscosity(i,  j+1,k+1)
                
                mu_plushalf_y_val  = (viscosity(i,j+1,k+1)+viscosity(i+1,j+1,k+1)+viscosity(i,j+2,k+1)+viscosity(i+1,j+2,k+1))/4.d0
                mu_minushalf_y_val = (viscosity(i,j+1,k+1)+viscosity(i+1,j+1,k+1)+viscosity(i,j  ,k+1)+viscosity(i+1,j  ,k+1))/4.d0
                
                mu_plushalf_z_val  = (viscosity(i,j+1,k+1)+viscosity(i+1,j+1,k+1)+viscosity(i,j+1,k+2)+viscosity(i+1,j+1,k+2))/4.d0
                mu_minushalf_z_val = (viscosity(i,j+1,k+1)+viscosity(i+1,j+1,k+1)+viscosity(i,j+1,k  )+viscosity(i+1,j+1,k  ))/4.d0
                
                rho_ugrid = (density(i+1,j+1,k+1)+density(i,j+1,k+1))/2.d0
                
                diag_val = 1.d0+(ax/2.d0)*(1.d0/rho_ugrid)*(mu_plushalf_x_val+mu_minushalf_x_val)+&
                     (ay/2.d0)*(1.d0/rho_ugrid)*(mu_plushalf_y_val+mu_minushalf_y_val)+&
                     (az/2.d0)*(1.d0/rho_ugrid)*(mu_plushalf_z_val+mu_minushalf_z_val)
                
                residual = (1.d0/diag_val) * ( &
                     (ax/2.d0)*(1.d0/rho_ugrid)*(mu_plushalf_x_val*u3(i+1,j,k)+mu_minushalf_x_val*u3(i-1,j,k)) + &
                     (ay/2.d0)*(1.d0/rho_ugrid)*(mu_plushalf_y_val*u3(i,j+1,k)+mu_minushalf_y_val*u3(i,j-1,k)) + &
                     (az/2.d0)*(1.d0/rho_ugrid)*(mu_plushalf_z_val*u_plusz    +mu_minushalf_z_val*u_minusz   ) + RHS(i,j,k))
                
                
                u3(i,j,k) = u3(i,j,k)+relax*(residual-u3(i,j,k))
                
             end do
          end do
       end do
       
    end if
    
    return
  end subroutine do_sor_u
  
  ! **********************************************************************************************
  
  subroutine do_sor_v(v3,RHS,viscosity,density,dx,dy,dz,dt,ex,ey,sx,sy,maxl,maxn,ex_max,flg,iteration_time)
    implicit none
    
    integer :: sx,sy,ex,ey,flg
    integer :: i,j,k,maxl,maxn,iteration_time,ex_max
    
    double precision :: v3(sx-1:ex+1,sy-1:ey+1,0:maxn-2), RHS(sx-1:ex+1,sy-1:ey+1,0:maxn-2)
    double precision :: viscosity(sx-2:ex+2,sy-2:ey+2,0:maxn)
    double precision :: density(sx-1:ex+1,sy-1:ey+1,0:maxn)
    double precision :: dx,dy,dz,dt,diag_val,relax,residual,ax,ay,az
    double precision :: mu_plushalf_x_val,mu_minushalf_x_val,mu_plushalf_y_val,mu_minushalf_y_val, &
         mu_plushalf_z_val,mu_minushalf_z_val
    
    double precision :: v_minusz,v_plusz,rho_vgrid
    
    relax=1.2d0
    ax=dt/(dx*dx)
    ay=dt/(dy*dy)
    az=dt/(dz*dz)
    
    if(mod(iteration_time,2).eq.0)then
       
       do k=0,maxn-2
          do j=sy,ey
             do i=mod(k+j+flg,2)+sx,ex,2
                
                if(k.eq.0)then
                   v_minusz=-2.d0*v3(i,j,0)+v3(i,j,1)/3.d0
                else
                   v_minusz=v3(i,j,k-1)
                end if
                
                if(k.eq.(maxn-2))then
                   v_plusz=-2.d0*v3(i,j,maxn-2)+v3(i,j,maxn-3)/3.d0
                else
                   v_plusz=v3(i,j,k+1)
                end if
                
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
                     (ax/2.d0)*(1.d0/rho_vgrid)*(mu_plushalf_x_val*v3(i+1,j,k)+mu_minushalf_x_val*v3(i-1,j,k)) + &
                     (ay/2.d0)*(1.d0/rho_vgrid)*(mu_plushalf_y_val*v3(i,j+1,k)+mu_minushalf_y_val*v3(i,j-1,k)) + &
                     (az/2.d0)*(1.d0/rho_vgrid)*(mu_plushalf_z_val*v_plusz    +mu_minushalf_z_val*v_minusz   ) + RHS(i,j,k))
                
                v3(i,j,k) = v3(i,j,k)+relax*(residual-v3(i,j,k))
                
             end do
          end do
       end do
       
    else
       
       do k=0,maxn-2
          do j=sy,ey
             do i=mod(j+flg,2)+sx,ex,2
                
                if(k.eq.0)then
                   v_minusz=-2.d0*v3(i,j,0)+v3(i,j,1)/3.d0
                else
                   v_minusz=v3(i,j,k-1)
                end if
              
                if(k.eq.(maxn-2))then
                   v_plusz=-2.d0*v3(i,j,maxn-2)+v3(i,j,maxn-3)/3.d0
                else
                   v_plusz=v3(i,j,k+1)
                end if
                
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
                     (ax/2.d0)*(1.d0/rho_vgrid)*(mu_plushalf_x_val*v3(i+1,j,k)+mu_minushalf_x_val*v3(i-1,j,k)) + &
                     (ay/2.d0)*(1.d0/rho_vgrid)*(mu_plushalf_y_val*v3(i,j+1,k)+mu_minushalf_y_val*v3(i,j-1,k)) + &
                     (az/2.d0)*(1.d0/rho_vgrid)*(mu_plushalf_z_val*v_plusz    +mu_minushalf_z_val*v_minusz   ) + RHS(i,j,k))
                
                v3(i,j,k) = v3(i,j,k)+relax*(residual-v3(i,j,k))
                
             end do
          end do
       end do
       
    end if
    
    return
  end subroutine do_sor_v
  
  ! **********************************************************************************************
  
  subroutine do_sor_w(w3,RHS,viscosity,density,dx,dy,dz,dt,ex,ey,sx,sy,maxl,maxn,ex_max,flg,iteration_time)
    implicit none
    
    integer :: sx,sy,ex,ey,flg
    integer :: i,j,k,maxl,maxn,iteration_time,ex_max
    
    double precision :: w3(sx-1:ex+1,sy-1:ey+1,0:maxn-1), RHS(sx-1:ex+1,sy-1:ey+1,0:maxn-1)
    double precision :: viscosity(sx-2:ex+2,sy-2:ey+2,0:maxn)
    double precision :: density(sx-1:ex+1,sy-1:ey+1,0:maxn)
    double precision :: dx,dy,dz,dt,diag_val,relax,residual,ax,ay,az, &
         mu_plushalf_x_val,mu_minushalf_x_val,mu_plushalf_y_val,mu_minushalf_y_val, &
         mu_plushalf_z_val,mu_minushalf_z_val
    
    double precision :: rho_wgrid
    
    relax=1.2d0
    ax=dt/(dx*dx)
    ay=dt/(dy*dy)
    az=dt/(dz*dz)
    
    if(mod(iteration_time,2).eq.0)then
       
       do k=1,maxn-2
          do j=sy,ey
             do i=mod(k+j+flg,2)+sx,ex,2
                
                mu_plushalf_x_val= (viscosity(i,j+1,k)+viscosity(i,j+1,k+1)+viscosity(i+1,j+1,k)+viscosity(i+1,j+1,k+1))/4.d0
                mu_minushalf_x_val=(viscosity(i,j+1,k)+viscosity(i,j+1,k+1)+viscosity(i-1,j+1,k)+viscosity(i-1,j+1,k+1))/4.d0
                
                mu_plushalf_y_val= (viscosity(i,j+1,k)+viscosity(i,j+1,k+1)+viscosity(i,j+2,k)+viscosity(i,j+2,k+1))/4.d0
                mu_minushalf_y_val=(viscosity(i,j+1,k)+viscosity(i,j+1,k+1)+viscosity(i,j  ,k)+viscosity(i,j  ,k+1))/4.d0
                
                mu_plushalf_z_val=  viscosity(i,j+1,k+1)
                mu_minushalf_z_val= viscosity(i,j+1,k  )
                
                rho_wgrid=(density(i,j+1,k+1)+density(i,j+1,k))/2.d0
                
                diag_val=1.d0+(ax/2.d0)*(1.d0/rho_wgrid)*(mu_plushalf_x_val+mu_minushalf_x_val)+&
                     (ay/2.d0)*(1.d0/rho_wgrid)*(mu_plushalf_y_val+mu_minushalf_y_val)+&
                     (az/2.d0)*(1.d0/rho_wgrid)*(mu_plushalf_z_val+mu_minushalf_z_val)
                
                residual = (1.d0/diag_val) * ( &
                     (ax/2.d0)*(1.d0/rho_wgrid)*(mu_plushalf_x_val*w3(i+1,j,k)+mu_minushalf_x_val*w3(i-1,j,k)) + &
                     (ay/2.d0)*(1.d0/rho_wgrid)*(mu_plushalf_y_val*w3(i,j+1,k)+mu_minushalf_y_val*w3(i,j-1,k)) + &
                     (az/2.d0)*(1.d0/rho_wgrid)*(mu_plushalf_z_val*w3(i,j,k+1)+mu_minushalf_z_val*w3(i,j,k-1)) + RHS(i,j,k))
                
                w3(i,j,k) = w3(i,j,k)+relax*(residual-w3(i,j,k))
             end do
          end do
       end do
       
    else
       
       do k=1,maxn-2
          do j=sy,ey
             do i=mod(j+flg,2)+sx,ex,2
                
                mu_plushalf_x_val= (viscosity(i,j+1,k)+viscosity(i,j+1,k+1)+viscosity(i+1,j+1,k)+viscosity(i+1,j+1,k+1))/4.d0
                mu_minushalf_x_val=(viscosity(i,j+1,k)+viscosity(i,j+1,k+1)+viscosity(i-1,j+1,k)+viscosity(i-1,j+1,k+1))/4.d0
                
                mu_plushalf_y_val= (viscosity(i,j+1,k)+viscosity(i,j+1,k+1)+viscosity(i,j+2,k)+viscosity(i,j+2,k+1))/4.d0
                mu_minushalf_y_val=(viscosity(i,j+1,k)+viscosity(i,j+1,k+1)+viscosity(i,j  ,k)+viscosity(i,j  ,k+1))/4.d0
                
                mu_plushalf_z_val=  viscosity(i,j+1,k+1)
                mu_minushalf_z_val= viscosity(i,j+1,k  )
                
                rho_wgrid=(density(i,j+1,k+1)+density(i,j+1,k))/2.d0
                
                diag_val=1.d0+(ax/2.d0)*(1.d0/rho_wgrid)*(mu_plushalf_x_val+mu_minushalf_x_val)+&
                     (ay/2.d0)*(1.d0/rho_wgrid)*(mu_plushalf_y_val+mu_minushalf_y_val)+&
                     (az/2.d0)*(1.d0/rho_wgrid)*(mu_plushalf_z_val+mu_minushalf_z_val)
                
                residual = (1.d0/diag_val) * ( &
                     (ax/2.d0)*(1.d0/rho_wgrid)*(mu_plushalf_x_val*w3(i+1,j,k)+mu_minushalf_x_val*w3(i-1,j,k)) + &
                     (ay/2.d0)*(1.d0/rho_wgrid)*(mu_plushalf_y_val*w3(i,j+1,k)+mu_minushalf_y_val*w3(i,j-1,k)) + &
                     (az/2.d0)*(1.d0/rho_wgrid)*(mu_plushalf_z_val*w3(i,j,k+1)+mu_minushalf_z_val*w3(i,j,k-1)) + RHS(i,j,k))
                
                w3(i,j,k) = w3(i,j,k)+relax*(residual-w3(i,j,k))
             end do
          end do
       end do
       
    end if
    
    return
  end subroutine do_sor_w
  

end module tpls_helmholtz
