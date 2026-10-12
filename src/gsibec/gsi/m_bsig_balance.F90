module m_bsig_balance
!$$$  module documentation block
!                .      .    .                                       .
! module:    m_bsig_balance
!
! abstract: variance added to t, vp and ps by the static B balance
!           (streamfunction -> t, vp, ps; unbalanced t -> ps when fut2ps).
!
!   bkerror applies B = Bal (S V H V S) Bal^T, where S holds dssv, V is the
!   vertical recursive filter (frfhvo) and H the horizontal one. The balance
!   acts column by column, so with the horizontal correlation diagonal taken
!   as one (as the normalization intends), the variance it adds to a balanced
!   variable is
!
!     var = || V (dssv_src * c) ||^2,
!
!   summed over levels, where c is the column of balance coefficients for that
!   output and dssv_src the standard deviation of the source control variable
!   (including the vertical normalization). The unbalanced variables are
!   independent of each other, so contributions from different sources add.
!   Velocity potential depends on streamfunction at the same level only, so
!   its contribution reduces to (bvz * bsig_sf)^2.
!
!   The variances are computed on first use and kept until
!   bsig_balance_final. The strong constraint (TLNMC) is not included.
!
! program history log:
!   2026-10-12  shlyaeva - initial code
!
!$$$ end documentation block

  use m_kinds, only: r_kind,i_kind
  implicit none

  private
  public :: bsig_balance_var3d
  public :: bsig_balance_var2d
  public :: bsig_balance_final

  real(r_kind),allocatable,save :: balvar_t(:,:,:)
  real(r_kind),allocatable,save :: balvar_vp(:,:,:)
  real(r_kind),allocatable,save :: balvar_ps(:,:)
  logical,save :: computed_ = .false.
  logical,save :: active_ = .false.

contains

!--------------------------------------------------------
  subroutine bsig_balance_var3d(name, k, var)
  ! Add the balanced variance of 3D control variable `name` at level k to
  ! var(lat2,lon2). Variables other than t and vp are left unchanged.
  character(len=*), intent(in)    :: name
  integer(i_kind),  intent(in)    :: k
  real(r_kind),     intent(inout) :: var(:,:)

  call compute_()
  if (.not. active_) return
  select case (trim(name))
  case ('t')
     var = var + balvar_t(:,:,k)
  case ('vp')
     var = var + balvar_vp(:,:,k)
  end select
  end subroutine bsig_balance_var3d

!--------------------------------------------------------
  subroutine bsig_balance_var2d(name, var)
  ! Add the balanced variance of 2D control variable `name` to var(lat2,lon2).
  ! Variables other than ps are left unchanged.
  character(len=*), intent(in)    :: name
  real(r_kind),     intent(inout) :: var(:,:)

  call compute_()
  if (.not. active_) return
  if (trim(name) == 'ps') var = var + balvar_ps
  end subroutine bsig_balance_var2d

!--------------------------------------------------------
  subroutine bsig_balance_final()
  if (allocated(balvar_t))  deallocate(balvar_t)
  if (allocated(balvar_vp)) deallocate(balvar_vp)
  if (allocated(balvar_ps)) deallocate(balvar_ps)
  computed_ = .false.
  active_ = .false.
  end subroutine bsig_balance_final

!--------------------------------------------------------
  subroutine compute_()
  use constants, only: zero
  use gridmod, only: lat2,lon2,nsig,regional
  use control_vectors, only: cvars3d,cvars2d
  use mpeu_util, only: getindex
  use berror, only: dssv,bsig3d,fpsproj,fut2ps
  use balmod, only: agvz,wgvz,bvz,pput
  use balmod, only: agvk,wgvk,bvk,agvk_lm,f1,fstat,ke_vp
  implicit none

  integer(i_kind) :: n_sf,n_t,n_vp,n_ps,i,j,k,l
  real(r_kind),allocatable :: u(:,:,:)

  if (computed_) return
  computed_ = .true.

  ! Same condition as dobal in bkerror
  n_sf = getindex(cvars3d,'sf')
  n_t  = getindex(cvars3d,'t')
  n_vp = getindex(cvars3d,'vp')
  n_ps = getindex(cvars2d,'ps')
  active_ = n_sf>0 .and. n_t>0 .and. n_vp>0 .and. n_ps>0
  if (.not. active_) return

  allocate(balvar_t(lat2,lon2,nsig),balvar_vp(lat2,lon2,nsig),balvar_ps(lat2,lon2))
  allocate(u(lat2,lon2,nsig))
  balvar_t = zero
  balvar_vp = zero
  balvar_ps = zero

  ! Velocity potential: same-level projection of streamfunction
  if (regional) then
     do k=1,ke_vp
        balvar_vp(:,:,k) = (bvk(:,:,k)*bsig3d(:,:,k,n_sf))**2
     end do
  else
     do k=1,nsig
        do j=1,lon2
           balvar_vp(:,j,k) = (bvz(:,k)*bsig3d(:,j,k,n_sf))**2
        end do
     end do
  end if

  ! Temperature: column projection of streamfunction, one output level at a time
  do k=1,nsig
     do l=1,nsig
        if (regional) then
           if (fstat) then
              u(:,:,l) = agvk_lm(k,l)*f1(:,:)*dssv(:,:,l,n_sf)
           else
              u(:,:,l) = agvk(:,:,k,l)*dssv(:,:,l,n_sf)
           end if
        else
           do j=1,lon2
              u(:,j,l) = agvz(:,k,l)*dssv(:,j,l,n_sf)
           end do
        end if
     end do
     call frfhvo(u,n_sf)
     balvar_t(:,:,k) = sum(u**2,dim=3)
  end do

  ! Surface pressure: column projection of streamfunction
  if (regional) then
     do l=1,nsig
        u(:,:,l) = wgvk(:,:,l)*dssv(:,:,l,n_sf)
     end do
  else
     do l=1,nsig
        do j=1,lon2
           u(:,j,l) = wgvz(:,l)*dssv(:,j,l,n_sf)
        end do
     end do
     if (.not. fpsproj) u(:,:,nsig) = zero
  end if
  call frfhvo(u,n_sf)
  balvar_ps = sum(u**2,dim=3)

  if (.not. regional) then
     if (fpsproj) then
        ! Projection of unbalanced temperature
        if (fut2ps) then
           do l=1,nsig
              do j=1,lon2
                 u(:,j,l) = pput(:,l)*dssv(:,j,l,n_t)
              end do
           end do
           call frfhvo(u,n_t)
           balvar_ps = balvar_ps + sum(u**2,dim=3)
        end if
     else
        ! Lowest level of unbalanced velocity potential replaces the top level of streamfunction
        do j=1,lon2
           do i=1,lat2
              balvar_ps(i,j) = balvar_ps(i,j) + (wgvz(i,nsig)*bsig3d(i,j,1,n_vp))**2
           end do
        end do
     end if
  end if

  deallocate(u)
  end subroutine compute_

end module m_bsig_balance
