#
#
#  This script manages the NOSE versions
#
#
#  Author:       Tomas Mancal, mancal@karlov.mff.cuni.cz
#  Last Change:  2008/07/01
#
#puts "Info: UberNose reads command line: $argv";

set first_arg [lindex $argv 0];

#
# First argument of the script on command line is interpreted as an option
# if it is any of the following
#
# --version      (prints info about the current version of NOSE)
# --list         (lists available versions)
# --reset        (resets settings to the default, which is the same as "--set-version latest")   
#
# If any of the above options are used, the script prerforms a corresponding action and exits
#
#

global curnosever;
global legacycode;

##############################################################
# Procedures
##############################################################
proc get_legacy { d } {

	set nlist [split $d .];
	set nr1 [lindex $nlist 0]; 
	set nr2 [lindex $nlist 1];
	if { $nr1 == 0 && $nr2 < 5 } {  
	
		return 1;
		
	} else {
	
		return 0;
		
	}
	
}

proc create_default_settings { } {

	global env;
	global exdir;
	global nhome;
	
	set home $env(HOME);
	set dirname ".nose";
	set dr [file join $home $dirname];
	set fl [file join $dr ubernose_settings.tcl];

	set f [open $fl w+];
	
	set defver [file tail [file dirname $exdir]];
	set leg [get_legacy $defver];
	puts $f "set curnosever $defver"
	puts $f "set legacycode $leg";	

	close $f;
}

proc read_settings { } {

	global curnosever;
	global legacycode;
	global env;
	global nhome;
	
	set home $env(HOME);
	
	set dirname ".nose";
	set dr [file join $home $dirname];
	set fl [file join $dr ubernose_settings.tcl];

	source $fl;
	
}

proc create_settings_dir { } {

	global env;
	global nhome;
	
	set home $env(HOME);
	set dirname ".nose";
	set dr [file join $home $dirname];
	
	if { [file exists $dr] == 1 && [file isdirectory $dr] == 1 } {
	
		
	} else {
	
		file mkdir $dr;
		create_default_settings;

	}   

}

proc update_settings { } {

	global curnosever;
	global legacycode;
	
	global env;
	global nhome;
	
	set home $env(HOME);
	set dirname ".nose";
	set dr [file join $home $dirname];
	set fl [file join $dr ubernose_settings.tcl];

	set f [open $fl w+]
	
	puts $f "set curnosever $curnosever"
	puts $f "set legacycode $legacycode";	

	close $f;
	
	if { $legacycode == 1 } {
	
		set pref [info script];
		set pref [file dirname $pref];
		puts "Info: !!!!! To run legacy version of NOSE, some enviromental variables need to be set "
		puts "Info: !!!!!                                                                           "
		puts "Info: !!!!! Source the following file [file join $pref nose-home $curnosever nose_setenv.sh]"
		puts "Info: !!!!!                                                                           "
		puts "Info: !!!!! You might need to source this file in every shell you use, or set it to be "
		puts "Info: !!!!! sourced at the shell opening --- this is the legacy NOSE --- sorry!!! "
		
	
	}

}

proc test_version { } {

	global legacycode;
	global curnosever;
	global env;
	global nhome;
	
	set home $nhome;#$env(HOME);
	
	if { $legacycode == 1 } {
		puts "Info    : Using legacy lnose script, version $curnosever ";
    	set ex [file join $home $curnosever bin nose1];
	} else {
		puts "Info    : Using the nose2 script, version $curnosever ";
    	set ex [file join $home $curnosever bin nose2];
	}

	if { [file exists $ex] == 0 } {
	
		puts "Info    : This version of NOSE is not present on this system";
		exit;
	}

}


##############################################################
#  Execution starts here
##############################################################
set chk no;
if { [string compare $first_arg "--check"] == 0} {

    puts "Running in test mode";
    set ex [file join .. .. src scripts  nose2];
   
    puts "Info   : executable script: $ex ";

    set chk yes;

    if { [llength $argv] == 0 } {
	
	set cmd [list $ex];
	
    } else {
	
	set cmd [list $ex];
	#set cmd [concat $cmd $argv];

    } 


} else {

	create_settings_dir;
	if { [catch { read_settings }] } {
		create_default_settings;
		read_settings;
	} 

#
# Set NOSE version to be used
#
if { [string compare $first_arg "--set-version"] == 0 } {

	global curnosever;
	global legacycode;
	
	create_settings_dir;
	read_settings;

	set log $curnosever;
	set curnosever [lindex $argv 1];
	
	# test new version	
	set legacycode [get_legacy [lindex $argv 1]];
	
	puts "Info: switching from NOSE version $log to $curnosever";
	test_version;
	update_settings;
	read_settings;
	puts "Info: update succesfull"
	exit;

}

#
# Get information about currently used version of NOSE
#
if { [string compare $first_arg "--version" ] == 0 } {

	read_settings;
	puts "Info: Curretly used version of NOSE is $curnosever";
	exit;

}

#
# Resets to default
#
if { [string compare $first_arg "--reset" ] == 0 } {

    create_settings_dir;
	read_settings;
	set log $curnosever;
    create_default_settings;
    read_settings;
	puts "Info: switching from NOSE version $log to $curnosever";
	test_version;
	puts "Info: update succesfull"
	exit;

}

#
# Lists avaibable versions 
# 
if { [string compare $first_arg "--list" ] == 0 } {

	global env;
	
	puts "Info: available versions are listed below"
	set dir [file join $nhome];
	foreach f [glob -nocomplain [file join $dir *]] {
		
		set vnr [file tail $f];
		if { [string match nose_* $vnr] == 0 } {
			set leg [get_legacy $vnr];
			if { $leg == 0 } {
		 
				puts "Info: version $vnr";
		
			} else {
		
				puts "Info: version $vnr (legacy version)";
			}
		}
	
	}
	
	puts "Info: The listed versions are not guarranteed to run correctly."
	exit;

}


#
# Generate a default working directory for the current version of NOSE
#
if  { [string compare $first_arg "--generate-dir" ] == 0 } {

	read_settings;
	
	puts "Info: generating working directory"
	set dir [file join $nhome];
	set fname [file join $dir $curnosever wdir];
	puts $fname;
	file copy -force $fname . ;
		
	puts "Info: directory named ./wdir has been created."
	exit;
 
}

 
 

if { $legacycode == 1 } { 
	puts "Info   : Using legacy lnose script, version $curnosever ";
    set ex [file join $nhome $curnosever bin nose1];
} else {
	puts "Info   : Using the nose2 script, version $curnosever ";
    set ex [file join $nhome $curnosever bin nose2];
} 
puts "Info   : executable script: $ex ";

if { [llength $argv] == 0 } {

	set cmd [list $ex];

} else {

	set cmd [list $ex];
	set cmd [concat $cmd]; # $argv];

} 

} 


#
# Do not use Expect package 
# 
global use_expect_package;
set use_expect_package 1

if { [string compare $first_arg "--no-expect" ] == 0 } {

	set use_expect_package 0	
	
#	puts "Info: Not using Expect package."
#	puts $argv;
	set argv [lrange $argv 1 [ llength $argv ] ]
#	puts [lrange $argv 1 [ llength $argv ] ];
#	puts $argv;

}



#if { [string compare $chk yes ] == 0 } {
#
#    set cmd [concat tclsh $cmd];
#
#} 

#package require Expect;

#spawn  $cmd;
#interact;

source $cmd

#set d [open "| $cmd " r+  ]
#while {[gets $d line ] >= 0 } {
#    puts $line;
#}
#close $d;

