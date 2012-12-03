#! /bin/csh -f

echo "LAMRANK is $LAMRANK"

xterm -e gdb $*



#if ("$LAMRANK" == "0") then
#  gdb $*
#else
#  $*
#endif

exit 0 


