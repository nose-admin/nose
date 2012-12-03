!
! Some tracing macros
!

character(len=32) :: tr_buff

#define ON_ENTER(procename,llevel)  write(tr_buff,'(a,a)') "Entering -> ", procname; call write_log_message(tr_buff,llevel)

#define ON_LEAVE(procename,llevel)  write(tr_buff,'(a,a)') "Leaving  -> ", procname; call write_log_message(tr_buff,llevel)
