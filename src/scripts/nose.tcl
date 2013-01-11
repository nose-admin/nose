#
#
#
#
#
#
#
#
#

set use_expect_package 0;

if { $use_expect_package == 1 } {
	package require Expect
}

package provide NOSE_conf 0.1

namespace eval NOSE_conf {

	variable parallel_execution;
	variable conf_parallel_options -np 1;
	variable recordList;
	variable parametricRecordList;
	variable parametricRecordCount;
	variable parametricRecordName;
	variable recordTouch;
	variable defaults;
	variable printDevice;
	variable printUnit;
	
	variable input_file;
	
	variable e_units;
	
	variable wn_to_cf 0.000188365156731;
	variable wl_to_cf 1883.6515673;

	variable matlist;
	
	variable to_source;
	variable smth2source 0;
	
	#
	###########################################################
	#  Procedures
	###########################################################
	#
	# Adds a record into the list
	#
	proc addRecord { rname values } {
    	upvar $values vals
    	variable recordList;

    	touchRecord $rname;

    	set recordList($rname) $values;
    
        #puts "$rname $values";
    
	}
	
	#
	# Adds a parametric record
	#
	proc addParametricRecord { rname pname values } {
		variable parametricRecordList;
		variable parametricRecordCount;
		variable parametricRecordName;
		variable parametricNames;
		
		if { ![info exists parametricRecordCount($rname)] } {
			set parametricRecordCount($rname) 0;
	    }
		incr parametricRecordCount($rname);
		
		lappend parametricNames($rname) $pname;
		set parametricRecordList($rname,$pname) $values;
		set parametricRecordName($rname,$parametricRecordCount($rname)) $pname;
	
		#puts "$rname $pname $values";
	}
	
	#
	# Adds a string record to the list
	#
	proc addStrRecord { rname val } {
    
    	set l [string length $val]
    	set ll [list v $l c $val];
    	addRecord $rname $ll

		#puts "$rname $val";

	}
	
	#
	# Touches a record
	#
	proc touchRecord { rname } {
    	variable defaults;
    	variable recordTouch;

		# if the record does not exist set the marker to 2
    	if { ![info exists recordTouch($rname)] } {
			set recordTouch($rname) 2;
    	}
    
    	# if we are setting defaults - reset the touch
    	if { $defaults == on } {
			set recordTouch($rname) 2; 
    	} else {
			set recordTouch($rname) [expr $recordTouch($rname) - 1];
    	}
    
    	if { $recordTouch($rname) == 0 } {
			warning "Record Already Exists!";
    	}
	}
	
	#
	#  Prints the record
	#
	proc putRecord { rname } {
    	variable recordList;
    	
 	    puts_extern $rname;
    
    	if { [info exists recordList($rname)] } {
			set ll $recordList($rname);
			foreach a $ll {
	    		print $a;
	    		puts_extern $a;
			}
    	} else {
			error "Record $rname does not exist"
    	}
    
	}
	
	proc getRecord { rname } {
	
    	variable recordList;
	
    	if { [info exists recordList($rname)] } {
			set ll $recordList($rname);
			foreach a $ll {
				set ven $a;
			}
    	} else {
			error "Record $rname does not exist"
    	}
    	
    	return $ven;
	
	}
	
	#
	#
	#
	proc putParametricRecordName { rname nr } {
		variable parametricRecordName;
		
		set name $parametricRecordName($rname,$nr);
	
		set l [list v [string length $name] c $name];

	    puts_extern $rname;
	
		foreach a $l {
	    	print $a;
	    	puts_extern $a;
		}

	
	}
	
	#
	#  Prints a parametric record
	#
	proc putParametricRecord { rname } {
		variable parametricRecordList;
		variable parametricNames;
		
		if { [info exists parametricNames($rname)] } {
		
			set nis_ll [list v [string length $rname] c $rname];
			
			#print $nis_ll;                                          # record name
			#puts_extern $nis_ll;                                          
			foreach a $nis_ll {
	    		print $a;
	    		puts_extern $a;
			}
			
			set ll $parametricNames($rname);
			set nr_pnames [llength $ll];
			set nis_ll [list s i $nr_pnames];
			
			#print $nis_ll;                                           # nr. of subrecords
			#puts_extern $nis_ll;
			foreach a $nis_ll {
	    		print $a;
	    		puts_extern $a;
			}
			
			foreach a $ll {
			
				set par $parametricRecordList($rname,$a);
				
				set par [join $par " "];
								
				set nis_ll [list v [string length $a] c $a];
				
				foreach a $nis_ll {
	    			print $a;
	    			puts_extern $a;
				}
				#print $nis_ll;                                       # subrecord name
				#puts_extern $nis_ll;
				
				set nis_ll [list v [string length $par] c $par];     
				
				foreach a $nis_ll {
	    			print $a;
	    			puts_extern $a;
				}
				#print $nis_ll;                                       # subrecord values
				#puts_extern $nis_ll;
				
			}
		
		} else {
			error "Record $rname does not exist";
		}
	
	} 
	
	#
	#
	#
	proc putInit { args } {
    	upvar $args ar
    	variable printDevice
    	variable printUnit

    	set printDevice 0

		# here should be some testing of the spawned process
	    if { ![info exists ar] } {
			set printDevice 0;
			set printUnit 0;
    	} else {
			set printDevice 1;
			set printUnit $ar; 
    	}

    	print "START_NIS_VER_1.0"

	}
	
	#
	#
	#
	proc putEnd { } {
    	variable printDevice
    	variable printUnit
    	global done

    	print "e"
    	if { $printDevice == 1 } {
			flush $printUnit;
						
			#while {[gets $printUnit line] >= 0 } {
			#	#puts ahoj;
			#	puts $line;
			#}
			
			#close $printUnit;
    	}
	}
	
	
	
	#
	#
	#
	proc puts_extern { str } {
	
		variable matlist;
		
		set genmat 1;
		
		if { $genmat } {
		
			lappend matlist $str;
		
		}
	
	}
	
	
	#
	#
	#
	proc print { char } {
    	variable printDevice;
    	variable printUnit;
    	global   use_expect_package;    	

#		if { $use_expect_package == 1 } {
#	    	if { $printDevice == 0 } {
#    			exp_send -- $char;
#    			exp_send "\n";
#				#puts $char;
#    		} else {
#				exp_send -- $char;
#    			exp_send "\n";
#				puts $printUnit $char;
#				#puts $char;
#			}
#    	} else {
	    	if { $printDevice == 0 } {
				puts $char;
    		} else {
				puts $printUnit $char;
				#puts $char;
			}
 #   	}
	}

	#
	#
	#
	proc checkReal { val comm nr } {


	}

	#
	#
	#
	proc checkInteger { val comm nr } {


	}
	
	#
	# Converts into the correct units
	#
	proc convert_energy { val } {
		variable e_units;
		variable wn_to_cf;
		variable wl_to_cf;
		
		switch -glob -- $e_units {
		
			
			wl { if { $val == 0.0 } { 
					return 0.0; 
				 } else {
				 	set ret [expr $wl_to_cf / $val ];
			        #puts  "wl: $val" 
				 	#puts  "cf: $ret" 
				 	return $ret
				 }                                      ;# wavelength } 
			wn { set ret [expr $wn_to_cf * $val ]
			     #puts  "wn: $val" 
				 #puts  "cf: $ret"
				 return $ret                            ;# wave number } 
			cf { 
			     #puts "cf: $val"
				 return $val                            ;# circular frequency 2*PI/fs } 
			default { error "Unknown energy units"; }
		
		} 	
	
	}
	
	 
	#
	########################################
	#  Input file commands
	########################################
	#
	#  moduleName
	#
	proc moduleName { mname } {
    	variable recordList; 

    	addStrRecord module_name $mname;

	}
	#
	#  moduleInpFile
	#
	proc moduleInpFile { mname } {
    	variable recordList; 

    	addStrRecord module_inp_file $mname;

	}
	#
	#  inputFile
	#
	proc inputFile { mname } {
    	variable recordList; 
		variable input_file;
		
		set input_file $mname;
				
    	addStrRecord inp_file $mname;

	}
	
	#
	#  gridStep
	#
	proc gridSteps { st1 st2 st3 } {  
    	variable recordList;

    	;# Check that they are reals
    	checkInteger $st1 gridStep 1
    	checkInteger $st2 gridStep 2
    	checkInteger $st3 gridStep 3

    	set ll [list v 3 i "$st1 $st2 $st3"];
    	addRecord grid_steps $ll;

	}
	#
	#  gridExtent
	#
	proc gridExtent { ge1 ge2 ge3 } {
    	variable recordList;

    	;# Check that they are integers
    	checkInteger $ge1 gridExtent 1
    	checkInteger $ge2 gridExtent 2
    	checkInteger $ge3 gridExtent 3

    	set ll [list v 3 i "$ge1 $ge2 $ge3"];
    	addRecord grid_extent $ll;

	}

	#
	# extendedGridExtent
	#
	proc extendedGridExtent { a } {
		variable recordList;
	
		checkInteger $a extendedGridExtent 1
		set ll [list m 1 1 i $a];
		addRecord extended_grid_extent $ll;
	
	}	

	#
	#  excPosition
	#
	proc excPosition { st1 st2 st3 st4 } {
		variable recordList;
	
		;# Check that they are integers
    	checkInteger $st1 excPosition 1
    	checkInteger $st2 excPosition 2
    	checkInteger $st3 excPosition 3
    	checkInteger $st4 excPosition 4
		
    	set ll [list v 4 i "$st1 $st2 $st3 $st4"];
		addRecord excPosition $ll;
	
	}

	#
	#  molecule 
	#
	proc molecule { a } {
		variable recordList;
	
		addStrRecord molecule $a;
	
	}

	#
	#  pathway
	#
	proc pathway { a } {
		variable recordList;
	
		addStrRecord pathway $a;
	
	}

	#
	#  runs 
	#
	proc runs { st1 st2} {
		variable recordList;
	
    	checkInteger $st1 runs 1;
    	checkInteger $st2 runs 2;

		set ll [list v 2 i "$st1 $st2"];
		addRecord runs $ll;
	
	}


	#
	#  spec_wini-dw-wst
	#
	proc spec_wini-dw-wst { st1 st2 st3 } {
		variable recordList;
	
		;# Check that they are integers
    	checkInteger $st1 spec_wini-dw-wst 1
    	checkInteger $st2 spec_wini-dw-wst 2
    	checkInteger $st3 spec_wini-dw-wst 3
		
    	set ll [list v 3 i "$st1 $st2 $st3"];
		addRecord spec_wini-dw-wst $ll;
	
	}
	#
	#  localBasis
	#
	proc localBasis { a } {
		variable recordList;
	
		addStrRecord localBasis $a;
	
	}

	#
	#  relaxation
	#
	proc relaxation { a } {
		variable recordList;
	
		addStrRecord relaxation $a;
	
	}


	#
	#  secular approximation
	#
	proc secular { a } {
		variable recordList;
	
		addStrRecord secular $a;
	
	}

	#
	#  dephasing
	#
	proc dephasing { a } {
		variable recordList;
	
		addStrRecord dephasing $a;
	
	}

	#
	#  feeding
	#
	proc feeding { a } {
		variable recordList;

		checkReal $a feeding 1
		set ll [list m 1 1 r $a];
		addRecord feeding $ll;
	
	}

	#
	#  draining
	#
	proc draining { a } {
		variable recordList;
	
		checkReal $a draining 1
		set ll [list m 1 1 r $a];
		addRecord draining $ll;
	
	}

	#
	#  parallel
	#
	proc parallel { a } {
		variable parallel_execution;
	
		set parallel_execution 0;
		if { $a == "yes" } {
			set parallel_execution 1;
		}
		addStrRecord parallel $a;
	}
	
	#
	# parallel options
	#
	proc parallel_options { a } {
	
		variable conf_parallel_options;
	
		set conf_parallel_options $a;

	}
	
	#
	#  tau of tau-dependent projector
	#
	proc tau { a } {
		variable recordList;
	
		checkReal $a tau 1
		set ll [list m 1 1 r $a];
		addRecord tau_projector $ll;
	
	}	
	#
	#  debug_gamma for montecarlo module
	#
	proc dgamma { a } {
		variable recordList;
	
		checkReal $a dgamma 1
		set ll [list m 1 1 r $a];
		addRecord debug_gamma $ll;
	
	}		
	#
	#  temperature
	#
	proc temperature { a } {
		variable recordList;
	
		checkReal $a temperature 1
		set ll [list m 1 1 r $a];
		addRecord temp $ll;

	}
	#
	# timeStep
	#
	proc timeStep { a } {
		variable recordList;
	
		checkReal $a timeStep 1
		set ll [list m 1 1 r $a];
		addRecord timeStep $ll;
	
	}
	#
	# logLevel
	#
	proc logLevel { a } {
		variable recordList;
	
		checkInteger $a logLevel 1
		set ll [list s i $a];
		addRecord logLevel $ll;
	
	}
	#
	#  pulsePolarization
	#
	proc pulsePolarization { ge1 ge2 ge3 } {
    	variable recordList;

    	;# Check that they are integers
    	checkReal $ge1 pulsePolarization 1
    	checkReal $ge2 pulsePolarization 2
    	checkReal $ge3 pulsePolarization 3

    	set ll [list v 3 r "$ge1 $ge2 $ge3"];
    	addRecord pulsePolarization $ll;

	}
	#
	# detectionPolarization
	#
	proc detectionPolarization { a } {
		variable recordList;
	
		checkReal $a detecttionPolarization 1
		set ll [list s r $a];
		addRecord detectionPolarization $ll;
	
	}


	# output
	
	#
	# outputName 
	#
	proc outputName { oname } {
    	variable recordList; 
		
    	addStrRecord output_name $oname;
	
	}
	
	#
	# outputDir
	#
	proc outputDir { mname } {
    	variable recordList; 
	
		if { ![file exists $mname] } {
			file mkdir $mname;
		} elseif { ![file isdirectory $mname] } {
			error "$mname is not a directory";
			exit;
		}
    	addStrRecord output_dir $mname;
	
	}
	
	#
	# outLogIds
	#
	proc outLogIds { ids } {
#    	variable recordList; 
#	
#		if { ![file exists $mname] } {
#			file mkdir $mname;
#		} elseif { ![file isdirectory $mname] } {
#			error "$mname is not a directory";
#			exit;
#		}
#   	addStrRecord output_dir $mname;
#	
	
	}


	#
	# output
	#
	proc output { args } {
		variable parametericRecordList;
		variable smth2source;
		variable to_source;
		
		set rstring "";
		set output_name [lindex $args 0];
		set len [llength $args];
		if { $len == 1 } {
			lappend args default;
			incr len;
		}
		set par [lrange $args 1 [expr $len - 1]];
		foreach a $par {
			lappend rstring $a;
		}
		
    	addParametricRecord output $output_name $rstring;
    	
    	if { $output_name == "dens_exc_block" } {
    
    		
    		if { $rstring == "default" } {
    			set smth2source 0;
    		} else {
    			set smth2source 1;
    		}
    		
    		set to_source $rstring;
    		
    	
    	} 
	
	}
	
	#
	# configuration of parametric record name
	#
	proc config { args } {
	
		
	
	}
	
	#
	# outLog
	#
	proc outLog { args } {
		variable parametericRecordList;
		
		set rstring "";
		set outlog_name [lindex $args 0];
		set len [llength $args];
		set par [lrange $args 1 [expr $len - 1]];
		foreach a $par {
			append rstring " " $a;
		}
		
    	addParametricRecord outLog $outlog_name $rstring;
	
    }

	#
	#
	#
	proc saveNIS { nis } {
		variable recordList;
		
		addStrRecord saveNIS $nis;
	
	}

	#
	#
	#
	proc saveGoft { nis } {
		variable recordList;
		
		addStrRecord saveGoft $nis;
	
	}

	#
	#
	#
	proc restartFreq { fr } {
		variable recordList;
		
		checkInteger $fr restartFreq 1;
		set rec [list s i $fr];
		
		addRecord restartFreq $rec; 
	
	}

	#
	#
	#
	proc realizations { av } {
		variable recordList;
		
		checkInteger $av realizations 1;
		set rec [list s i $av];
		
		addRecord realizations $rec; 
	
	}

	#
	# Sets the NOSE temporary directory
	#
	proc tempDir { td } {
		variable recordList;
		
		addStrRecord tempDir $td;
	
	}

	#
	# Determines if the temporary directory should be
	# deleted after the run.
	#
	proc deleteTemp { td } {
		variable recordList;
		
		addStrRecord deleteTemp $td;
	
	}

	#
	# Sets current units
	#
	proc units { args } {
		variable e_units;
		
		set e_units [lindex $args 0];
	
	}
	
	#
	# Pulse frequencies 
	#
	proc pulseFrequency { args } {
		variable e_units;
		
		# convert to correct units
		
		set f1 [convert_energy 0.0];
		set f2 [convert_energy 0.0];
		set f3 [convert_energy 0.0];
		
		checkReal $f1 pulseFrequency 1
		checkReal $f2 pulseFrequency 2
		checkReal $f3 pulseFrequency 3
		
		set rec [list v 3 r $f1 $f2 $f3]; 
		
		addRecord pulseFrequencies $rec;
	
	}
	
	#
	# RWA frequency
	#
	proc rwa { a } {
		variable recordList;
		variable e_units;
		
		checkReal $a rwa 1;
		
		set ff [convert_energy $a];
		
		set ll [list s r $ff];
		
		addRecord rwa $ll;
	
	}

	#
	# Module method
	#
	proc moduleMethod { m } {
	
		variable recordList;
		
		addStrRecord moduleMethod $m;
	
	}

	#
	# completeResultFile 
	#
	proc completeResultFile { oname } {
    	variable recordList; 
		
	    addStrRecord completeResultFile $oname;
	
	}

	#
	##############################################
	#  End Input Commands
	##############################################


	##############################################
	#  Set defauls
	##############################################
	set defaults on;

	inputFile        nose.ssf
	moduleInpFile    none   ;# !!! what to do with "/" chacter?
	gridSteps        1 100 1
	gridExtent       1000 1 1000
	parallel         no
	logLevel         5
	saveGoft         no
	saveNIS          no
	outputDir        ./
	realizations     1
	restartFreq      100
	parallel         no
	parallel_options "-np 1"
	
	moduleMethod PT2-RC
	excPosition       1 4 2 2 
	molecule         dimer
	pathway          R1
	runs             1  -50
	spec_wini-dw-wst 12100 100 1          
	localBasis       no
	relaxation       no
	secular          yes
	dephasing        yes
	feeding           0.
	draining          0.
    
	set defaults off;
	##############################################
	#  End setting defaults
	##############################################

	namespace export moduleName;
	namespace export moduleInpFile;
	namespace export inputFile
	namespace export gridSteps;
	namespace export gridExtent;
	namespace export extendedGridExtent;
	namespace export parallel;
	namespace export tau;
	namespace export dgamma;
	namespace export temperature;
	namespace export timeStep;
	namespace export logLevel;
	namespace export detectionPolarization;
	namespace export pulsePolarization;
	namespace export outputName;
	namespace export outputDir;
	namespace export output;
	namespace export outLog;
	namespace export outLogIds;
	namespace export saveNIS;
	namespace export saveGoft;
	namespace export restartFreq;
	namespace export realizations;
	namespace export tempDir;
	namespace export deleteTemp;
	namespace export units;
	namespace export pulseFrequency;
	namespace export putParametricRecordName;
	namespace export rwa;
	namespace export moduleMethod;
	namespace export completeResultFile;
	namespace export parallel_options;
	namespace export excPosition;
	namespace export molecule;
	namespace export pathway;
	namespace export runs;
	namespace export spec_wini-dw-wst;
	namespace export localBasis;
	namespace export relaxation;
	namespace export secular;
	namespace export dephasing;
	namespace export feeding;
	namespace export draining;

} ;# END namespace NOSE_conf


package provide NOSE_ssf 0.1

namespace eval NOSE_ssf {

	#
	# Syntax checking
	#
	variable level 0
	variable level_name [list GROUND];
	
	variable block_indices [list];
	variable current_block_id 0;
	variable block_count 0;
	variable goft_count  0;

	variable transition_ids;
	variable transition_coordinates;
	variable transition_dipoles;
	variable transition_dlengths;
	variable transition_energies;
	variable transition_gofts;
    variable transition_count;
    variable transition_disorders;
	variable couplings;
	variable int_couplings;
	variable current_block_id_1 0;
	variable current_block_id_2 0;
	variable mode_id 0;
	variable goft_id 0;
	variable mode_type;
	variable mode_count 0;
	variable reorganization_energy;
    variable qo_levels;
    variable qo_frequency;
    variable qo_reorganization_energy;
    variable qo_huang_rhys_factor;
    variable correlation_time;
    variable rstrength;
    
	variable mcounts;
	variable nose_home;
	variable nexdir;
	
	variable module_method;
	variable periodicity_dim 3;
	variable periodicity;
	variable lattice_cell;
	variable cell_coupling;
	variable lattice_vec_r1;
	variable lattice_vec_r2;
	variable lattice_vec_r3;
	variable periodic 0;
    
    variable basis_type
    variable current_ic_id;
    variable dm_nr_states;
    variable dm_values_r;
    variable dm_values_i;
    variable ic_count 0;
    
    set nose_home $nhome ;#$env(NOSE_HOME);
    set nexdir $exdir;
    
	proc check_level { lev lname } {
		variable level;
		variable level_name;
		
		set last_level [lindex $level_name [expr [llength $level_name] - 1]];
		if { $lev != $level } {
			return -1;
		} 
		
		if { ![string match $lname $last_level] } {
			return 1;
		} 
	
		return 0;
		
	}

	proc print_LEVEL { } {
		variable level;
		variable level_name;
	
		puts "level: $level - $level_name";
	
	}

	proc print_ERROR { who } {
	
		puts "Input file -> Syntax error: $who";
		error $who;
	
	}
	
	#
	# Increase level
	#
	proc incr_level { lname } {
		variable level;
		variable level_name;
	
		incr level;
		lappend level_name $lname;
	
	}
	
	#
	# Decrease level
	#
	proc decr_level { } {
		variable level;
		variable level_name;

		set le [llength $level_name];
		set level_name [lrange $level_name 0 [expr $le - 2]];
		set level [expr $level - 1];
	
	}




    #
	#  Input file commands
    #

	proc BEGIN_SSF { args } {
		variable level;
		variable level_name;
		variable lattice_cell;			
		variable lattice_vec_r1;
		variable lattice_vec_r2;
		variable lattice_vec_r3;
		
		
		if {[check_level 0 GROUND] == 0 } {
		
			incr_level SSF;		
			
			set lattice_cell(1) 0;
			set lattice_cell(2) 0;
			set lattice_cell(3) 0;
			set lattice_vec_r1(1) 0;
			set lattice_vec_r1(2) 0;
			set lattice_vec_r1(3) 0;
			set lattice_vec_r2(1) 0;
			set lattice_vec_r2(2) 0;
			set lattice_vec_r2(3) 0;
			set lattice_vec_r3(1) 0;
			set lattice_vec_r3(2) 0;
			set lattice_vec_r3(3) 0;
		
		
		} else {
		
			print_ERROR BEGIN_SSF;
			exit;
		
		}

	}

	proc END_SSF { } {
		variable level;
		variable level_name;

		if {[check_level 1 SSF] == 0 } {

			decr_level;			
						
		} else {
			
			print_ERROR END_SSF;
			exit;
		
		}
	}

	proc BEGIN_BLOCK { block_id } {
		variable level;
		variable level_name;

		variable block_count;
		variable block_indices;
		variable current_block_id;
		variable transition_ids;
		variable transition_count;

		if {[check_level 1 SSF] == 0 } {
		
			incr_level BLOCK;		

			#
			# test block id
			#
			
			if { $block_id <= 0 } {
			
				print_ERROR "Block id must be positive!";
				exit;
			
			}
			
			set ind [lsearch $block_indices $block_id];			
			if { $ind >= 0 } {
			
				print_ERROR "Block id $block_id already exists";
				exit;
			
			}

			# save block id
			
			lappend block_indices $block_id
			set current_block_id  $block_id
			set transition_ids($current_block_id) [list]; 
			set transition_count($current_block_id) 0;
			incr block_count;
		
		} else {
		
			print_ERROR BEGIN_BLOCK;
			exit;
		
		}

	}

	proc END_BLOCK { } {
		variable level;
		variable level_name;
		variable current_block_id

		if {[check_level 2 BLOCK] == 0 } {

			decr_level;			
						
			set current_block_id  0


		} else {
			
			print_ERROR END_BLOCK;
			exit;
		
		}

	}

    #
	#  Specifies periodicity of the problem
    #
	proc BEGIN_PERIODIC { args } {
		variable level;
		variable level_name;

		variable periodic;

		if {[check_level 2 BLOCK] == 0 } {
		
			incr_level PERIODIC;		
			set periodic 1;
		
		} else {
		
			print_ERROR BEGIN_PERIODIC;
			exit;
		
		}

	}


	proc END_PERIODIC { } {
		variable level;
		variable level_name;
		variable lattice_cell;

		if {[check_level 3 PERIODIC] == 0 } {

			decr_level;			
			set lattice_cell(1) 0;
			set lattice_cell(2) 0;
			set lattice_cell(3) 0;				
			
		} else {
			
			print_ERROR END_PERIODIC;
			exit;
		
		}

	}



	#
	# Optical transition
	#
	proc TRANSITION { args } {
		variable level;
		variable level_name;

		variable current_block_id;
		variable transition_coordinates;
		variable transition_dipoles;
		variable transition_dlengths;
		variable transition_ids;
		variable transition_energies;
		variable transition_gofts;
		variable transition_disorders;
		variable transition_vibrational_levels;
		variable transition_count;

		if {[check_level 2 BLOCK] == 0 } {
		
			# check number of arguments
			set len [llength $args];
			if { $len < 9 } {
			
				print_ERROR "Wrong number of arguments in TRANSITION";
				exit;
			
			}
			
			set id [lindex $args 0]

			set tr_indices $transition_ids($current_block_id);
			
			set ind [lsearch $tr_indices $id];			
			if { $ind >= 0 } {
			
				print_ERROR "Transition id $id already exists in block $block_id";
				exit;
			
			}
			
			lappend transition_ids($current_block_id) $id;
						
			set transition_coordinates($current_block_id,$id,1) [lindex $args 1];
			set transition_coordinates($current_block_id,$id,2) [lindex $args 2];
			set transition_coordinates($current_block_id,$id,3) [lindex $args 3];
			set transition_dipoles($current_block_id,$id,1)     [lindex $args 4];
			set transition_dipoles($current_block_id,$id,2)     [lindex $args 5];
			set transition_dipoles($current_block_id,$id,3)     [lindex $args 6];
			set transition_dlengths($current_block_id,$id)      [lindex $args 7];
			set transition_energies($current_block_id,$id)      [NOSE_conf::convert_energy  [lindex $args 8]];
			set transition_gofts($current_block_id,$id)         [lindex $args 9];

			if { $len > 10 } {
			
				set transition_disorders($current_block_id,$id,$id)         [lindex $args 10];				
			
			} else {
			
				set transition_disorders($current_block_id,$id,$id)         0;
				
			}
			
			if { $len > 11 } {
						
				set transition_vibrational_levels($current_block_id,$id)         [lindex $args 11];				
						
			} else {
						
				set transition_vibrational_levels($current_block_id,$id)         0;
					
			}
			
			incr transition_count($current_block_id);
			
		} else {
			
			print_ERROR TRANSITION;
			exit;
		
		}

	}


	#
	# Dimension of the periodicity
	#
	proc DIMENSIONS { args } {
		variable level;
		variable level_name;
		
		variable periodicity_dim;

		if {[check_level 3 PERIODIC] == 0 } {
	
			# check number of arguments
			set len [llength $args];
			if { $len != 1 } {
			
				print_ERROR "Wrong number of arguments in DIMENSIONS";
				exit;
			
			}
			
			set periodicity_dim [lindex $args 0];
			
#			puts "$periodicity_dim";
		
		
		} else {
			
			print_ERROR DIMENSIONS;
			exit;
		
		}
		
	}

	#
	# Length of the periodicity
	#
	proc PERIODICITY { args } {
		variable level;
		variable level_name;
		
		variable periodicity;

		if {[check_level 3 PERIODIC] == 0 } {
	
			# check number of arguments
			set len [llength $args];
			if { $len != 1 } {
			
				print_ERROR "Wrong number of arguments in PERIODICITY";
				exit;
			
			}
			
			set periodicity [lindex $args 0];
			
#			puts "$periodicity";
		
		} else {
			
			print_ERROR PERIODICITY;
			exit;
		
		}
		
	}

	#
	# Lattice vector 1
	#
	proc R1 { args } {
		variable level;
		variable level_name;
		
		variable lattice_vec_r1;
		variable periodicity_dim;

		if {[check_level 3 PERIODIC] == 0 } {
	
			# check number of arguments
			set len [llength $args];
			if { $len != 3 } {
			
				print_ERROR "Wrong number of arguments in R1";
				exit;
			
			}
			
			set lattice_vec_r1(1) [lindex $args 0];
			set lattice_vec_r1(2) [lindex $args 1];
			set lattice_vec_r1(3) [lindex $args 2];
			
		
		} else {
			
			print_ERROR R1;
		
		}
		
	}

	#
	# Lattice vector 2
	#
	proc R2 { args } {
		variable level;
		variable level_name;
		
		variable lattice_vec_r2;
		variable periodicity_dim;

		if {[check_level 3 PERIODIC] == 0 } {
	
			# check number of arguments
			set len [llength $args];
			if { $len != 3 } {
			
				print_ERROR "Wrong number of arguments in R2";
				exit;
			
			}
			
			set lattice_vec_r2(1) [lindex $args 0];
			set lattice_vec_r2(2) [lindex $args 1];
			set lattice_vec_r2(3) [lindex $args 2];
			
		
		} else {
			
			print_ERROR R2;
			exit;
		
		}
		
	}

	#
	# Lattice vector 3
	#
	proc R3 { args } {
		variable level;
		variable level_name;
		
		variable lattice_vec_r3;
		variable periodicity_dim;

		if {[check_level 3 PERIODIC] == 0 } {
	
			# check number of arguments
			set len [llength $args];
			if { $len != 3 } {
			
				print_ERROR "Wrong number of arguments in R3";
				exit;
			
			}
			
			set lattice_vec_r3(1) [lindex $args 0];
			set lattice_vec_r3(2) [lindex $args 1];
			set lattice_vec_r3(3) [lindex $args 2];
			


		} else {
			
			print_ERROR R3;
			exit;
		
		}
		
	}

	proc TURN { args } {
		variable level;
		variable level_name;
		
		variable lattice_cell;
		variable lattice_cell_turn;
		variable periodicity_dim;

		if {[check_level 3 PERIODIC] == 0 } {
	
			# check number of arguments
			set len [llength $args];
			if { $len != 3 } {
			
				print_ERROR "Wrong number of arguments in TURN";
				exit;
			
			}
			
			set lattice_cell_turn(1) [lindex $args 0];
			set lattice_cell_turn(2) 0;
			set lattice_cell_turn(3) 0;
			
			if {$periodicity_dim > 1} {
			
				set lattice_cell_turn(2) [lindex $args 1];
				
				if { $periodicity_dim > 2} {
				
					set lattice_cell_turn(3) [lindex $args 2];
				
				}
			}


		} else {
			
			print_ERROR TURN;
			exit;
		
		}
	
	
	
	
	}

	#
	# Descriptor of the cell to interact with
	#
	proc CELL { args } {
		variable level;
		variable level_name;
		
		variable lattice_cell;
		variable periodicity_dim;

		if {[check_level 3 PERIODIC] == 0 } {
	
			# check number of arguments
			set len [llength $args];
			if { $len != $periodicity_dim } {
			
				print_ERROR "Wrong number of arguments in CELL";
				exit;
			
			}
			
			set lattice_cell(1) [lindex $args 0];
			set lattice_cell(2) 0;
			set lattice_cell(3) 0;
			
			if {$periodicity_dim > 1} {
			
				set lattice_cell(2) [lindex $args 1];
				
				if { $periodicity_dim > 2} {
				
					set lattice_cell(3) [lindex $args 2];
				
				}
			}


		} else {
			
			print_ERROR CELL;
			exit;
		
		}
		
	}


	#
	# Excitonic coupling
	#
	proc COUPLING { args } {
		variable level;
		variable level_name;
		
		variable couplings;
		variable int_couplings;
		variable current_block_id;
		variable current_block_id_1;
		variable current_block_id_2;
		variable transition_disorders;
		variable lattice_cell;
		variable cell_coupling;

		if {[check_level 2 BLOCK] == 0 } {
		
			# check number of arguments
			set len [llength $args];
			
			if { $len == 3} {
			
				set transition_disorders($current_block_id,[lindex $args 0],[lindex $args 1]) 0;
				set transition_disorders($current_block_id,[lindex $args 1],[lindex $args 0]) 0;
    			set couplings($current_block_id,[lindex $args 0],[lindex $args 1]) [NOSE_conf::convert_energy [lindex $args 2]];
	    		set couplings($current_block_id,[lindex $args 1],[lindex $args 0]) [NOSE_conf::convert_energy [lindex $args 2]];
						
			} elseif { $len == 4} {
			
				set transition_disorders($current_block_id,[lindex $args 0],[lindex $args 1]) [lindex $args 3];
				set transition_disorders($current_block_id,[lindex $args 1],[lindex $args 0]) [lindex $args 3];
    			set couplings($current_block_id,[lindex $args 0],[lindex $args 1]) [NOSE_conf::convert_energy [lindex $args 2]];
	    		set couplings($current_block_id,[lindex $args 1],[lindex $args 0]) [NOSE_conf::convert_energy [lindex $args 2]];
			
			} else {
			
				print_ERROR "Wrong number of arguments in COUPLING";
				exit;			
			
			}
		
		} elseif {[check_level 2 INTER_BLOCK] == 0} {
		
			# check number of arguments
			set len [llength $args];
			if { $len != 3 } {
			
				print_ERROR "Wrong number of arguments in COUPLING";
				exit;
			
			}
			
			
			set int_couplings($current_block_id_1,$current_block_id_2,[lindex $args 0],[lindex $args 1]) [NOSE_conf::convert_energy [lindex $args 2]];
			
		} elseif { [check_level 3 PERIODIC] == 0} {
		
			# check number of arguments
			set len [llength $args];
			if { $len != 3 } {
			
				print_ERROR "Wrong number of arguments in COUPLING";
				exit;
				
			} else {
			
				set cell_coupling($lattice_cell(1),$lattice_cell(2),[lindex $args 0],[lindex $args 1]) [NOSE_conf::convert_energy [lindex $args 2]];				
			
			}
		
		} else {
			
			print_ERROR COUPLING;
			exit;
		
		}
	}

	#
	# Rotational strength
	#
	proc ROTATIONAL_STRENGTH { args } {
		variable level;
		variable level_name;
		variable rstrength;
		
		variable current_block_id;

		if {[check_level 2 BLOCK] == 0 } {
		
			# check number of arguments
			set len [llength $args];
			if { $len != 2 } {
			
				print_ERROR "Wrong number of arguments in ROTATIONAL_STRENGTH";
				exit;
			
			}
			
			set rstrength($current_block_id,[lindex $args 0]) [lindex $args 1];
			
		
		} else {
			
			print_ERROR ROTATIONAL_STRENGTH;
			exit;
		
		}
	}


	proc UNITS { args } {
		variable level;
		variable level_name;

		if { $level > 0 } {


		} else {
			
			print_ERROR UNITS;
			exit;
		
		}

	}

	proc UNITS_ENERGY { un } {
		variable level;
		variable level_name;

		if { $level > 0 } {

			set NOSE_conf::e_units $un;

		} else {
			
			print_ERROR UNITS_ENERGY;
			exit;
		
		}

	}

	proc UNITS_TIME { un } {
		variable level;
		variable level_name;

		if { $level > 0 } {


		} else {
			
			print_ERROR UNITS_TIME;
			exit;
		
		}

	}
	
	proc BEGIN_QUANTUM_OSCILLATOR { id } {
		variable level;
		variable level_name;
		variable qo_id;
		
		if {[check_level 1 SSF] == 0 } {
		
			incr_level QO_OSCILLATOR;		

			set qo_id $id;
		
		
		} else {
		
			print_ERROR BEGIN_QUANTUM_OSCILLATOR;
			exit;
		
		}

	}

	proc END_QUANTUM_OSCILLATOR { } {
		variable level;
		variable level_name;

		if {[check_level 2 QO_OSCILLATOR] == 0 } {

			decr_level;			
						


		} else {
			
			print_ERROR END_QUANTUM_OSCILLATOR;
			exit;
		
		}

	}

	proc QO_HUANG_RHYS_FACTOR { type } {
		variable level;
		variable level_name;
		variable qo_id;
		variable qo_huang_rhys_factor;
		variable qo_reorganization_energy;
		
		if { [info exists qo_reorganization_energy($qo_id) ] } {
			puts "cannot input both Huang-Rhys factor and reorganization energy";
			print_ERROR QO_HUANG_RHYS_FACTOR;
			exit;
		}

		if {[check_level 2 QO_OSCILLATOR] == 0 } {
		
			set qo_huang_rhys_factor($qo_id) $type 

		} else {
			
			print_ERROR QO_HUANG_RHYS_FACTOR;
			exit;
		
		}

	}
	
	proc QO_REORGANIZATION_ENERGY { type } {
		variable level;
		variable level_name;
		variable qo_id;
		variable qo_reorganization_energy;
		variable qo_huang_rhys_factor;
		
		if { [info exists qo_huang_rhys_factor($qo_id) ] } {
			puts "cannot input both Huang-Rhys factor and reorganization energy";
			print_ERROR QO_REORGANIZATION_ENERGY;
			exit;
		}
		

		if {[check_level 2 QO_OSCILLATOR] == 0 } {
		
			set qo_reorganization_energy($qo_id) [NOSE_conf::convert_energy $type]; 

		} else {
			
			print_ERROR QO_REORGANIZATION_ENERGY;
			exit;
		
		}

	}
	
	proc QO_LEVELS { type } {
		variable level;
		variable level_name;
		variable qo_id;
		variable qo_levels;

		if {[check_level 2 QO_OSCILLATOR] == 0 } {
		
			set qo_levels($qo_id) $type 

		} else {
			
			print_ERROR QO_LEVELS;
			exit;
		
		}

	}
	
	proc QO_FREQUENCY { type } {
		variable level;
		variable level_name;
		variable qo_id;
		variable qo_frequency;

		if {[check_level 2 QO_OSCILLATOR] == 0 } {
		
			set qo_frequency($qo_id) [NOSE_conf::convert_energy $type]; 

		} else {
			
			print_ERROR QO_FREQUENCY;
			exit;
		
		}

	}	
	
	proc QO_OSCILLATOR { qo_oscillator } {
		variable level;
		variable level_name;

		if {[check_level 3 GROUP] == 0 } {
		

		} else {
			
			print_ERROR QO_OSCILLATOR;
			exit;
		
		}

	} 	

	proc BEGIN_CORRF { id } {
		variable level;
		variable level_name;
		variable goft_count;
		variable goft_id;

		if {[check_level 1 SSF] == 0 } {
		
			incr_level CORRF;		

			# save corrf id
			set goft_id $id;
				
			incr goft_count;
		
		
		} else {
		
			print_ERROR BEGIN_CORRF;
			exit;
		
		}

	}
	
	proc BEGIN_MODE { id } {
		variable level;
		variable level_name;
		variable mode_count;
		variable mode_id;
		
		if {[check_level 2 CORRF] == 0 } {
			incr_level MODE;
			incr mode_count;		
			set mode_id $id;
			set mids($mode_count) $id;
		
		} else {
			print_ERROR BEGIN_MODE;
			exit;
		}

	}

	proc END_CORRF { } {
		variable level;
		variable level_name;
		variable goft_id;
		variable mode_count;
		variable mcounts;

		if {[check_level 2 CORRF] == 0 } {

			decr_level;			
			set mcounts($goft_id) $mode_count;					
		    set goft_id 0;
			set mode_count 0;
	
	
		} else {
			
			print_ERROR END_CORRF;
			exit;
		
		}

	}

	proc END_MODE { } {
		variable level;
		variable level_name;
		variable mode_id;

		if {[check_level 3 MODE] == 0 } {

			decr_level;			
			set mode_id 0;


		} else {
			
			print_ERROR END_MODE;
			exit;
		
		}
		
	}

	proc MODE_TYPE { type } {
		variable level;
		variable level_name;
		variable mode_type;
		variable goft_id;
		variable mode_id;

		if {[check_level 3 MODE] == 0 } {
		
			set mode_type($goft_id,$mode_id) $type

		} else {
			
			print_ERROR MODE_TYPE;
			exit;
		
		}
	}

	#
	#
	#
	proc REORGANIZATION_ENERGY { lambda } {
		variable level;
		variable level_name;
		variable reorganization_energy;
		variable goft_id;
		variable mode_id;

		if {[check_level 3 MODE] == 0 } {
		
			set reorganization_energy($goft_id,$mode_id) [NOSE_conf::convert_energy $lambda ];

		} else {
			
			print_ERROR REORGANIZATION_ENERGY;
			exit;
		
		}

	}

	#
	#
	#
	proc CORRELATION_TIME { tc } {
		variable level;
		variable level_name;
    	variable correlation_time;
		variable goft_id;
		variable mode_id;

		if {[check_level 3 MODE] == 0 } {
		
			set correlation_time($goft_id,$mode_id) $tc;

		} else {
			
			print_ERROR CORRELATION_TIME;
			exit;
		
		}

	}
	
	#
	#
	#
	proc BROWNIAN_OMEGA { tc } {
		variable level;
		variable level_name;
    	variable brownian_omega;
		variable goft_id;
		variable mode_id;

		if {[check_level 3 MODE] == 0 } {
		
			set brownian_omega($goft_id,$mode_id) [NOSE_conf::convert_energy $tc ];

		} else {
			
			print_ERROR BROWNIAN_OMEGA;
			exit;
		
		}

	}
	
	
	#
	#  Homogeneous width
	#
	proc GAMMA { tc } {
		variable level;
		variable level_name;
    	variable gamma;
		variable goft_id;
		variable mode_id;

		if {[check_level 3 MODE] == 0 } {
		
			set gamma($goft_id,$mode_id) $tc;

		} else {
			
			print_ERROR GAMMA;
			exit;
		
		}

	}

	proc BEGIN_DISORDER { id } {
		variable level;
		variable level_name;
		variable dis_id;

		if {[check_level 1 SSF] == 0 } {
		
			incr_level DISORDER;		

			set dis_id $id;
		
		
		} else {
		
			print_ERROR BEGIN_DISORDER;
			exit;
		
		}

	}

	proc END_DISORDER { } {
		variable level;
		variable level_name;

		if {[check_level 2 DISORDER] == 0 } {

			decr_level;			
						


		} else {
			
			print_ERROR END_DISORDER;
			exit;
		
		}

	}

	proc TYPE { type } {
		variable level;
		variable level_name;
		variable dis_type;
		variable dis_id;

		if {[check_level 2 DISORDER] == 0 } {
		
			set dis_type($dis_id) $type 

		} else {
			
			print_ERROR TYPE;
			exit;
		
		}

	}

	proc TARGETS { args } {
		variable level;
		variable level_name;

		if {[check_level 2 DISORDER] == 0 } {
		

		} else {
			
			print_ERROR TARGETS;
			exit;
		
		}

	}

	proc WIDTH { width } {
		variable level;
		variable level_name;
		variable dis_id;
		variable dis_width;

		if {[check_level 2 DISORDER] == 0 } {
		
			set dis_width($dis_id) [NOSE_conf::convert_energy $width] ;

		} else {
			
			print_ERROR WIDTH;
			exit;
		
		}

	}


	proc BEGIN_GROUP { id } {
		variable level;
		variable level_name;

		if {[check_level 2 BLOCK] == 0 } {
		
			incr_level GROUP;		

			# save group id
		
		
		} else {
		
			print_ERROR BEGIN_GROUP;
			exit;
		
		}

	} 

	proc END_GROUP { } {
		variable level;
		variable level_name;

		if {[check_level 3 GROUP] == 0 } {

			decr_level;			
						


		} else {
			
			print_ERROR END_GROUP;
			exit;
		
		}

	}

	proc MEMBERS { args } {
		variable level;
		variable level_name;

		if {[check_level 3 GROUP] == 0 } {
		

		} else {
			
			print_ERROR MEMBERS;
			exit;
		
		}

	}

	proc DISORDER { disorder } {
		variable level;
		variable level_name;

		if {[check_level 3 GROUP] == 0 } {
		

		} else {
			
			print_ERROR DISORDER;
			exit;
		
		}

	} 


	proc BEGIN_INTER_BLOCK { b1 b2 } {
		variable level;
		variable level_name;
		variable current_block_id_1;
		variable current_block_id_2;

		if {[check_level 1 SSF] == 0 } {
		
			incr_level INTER_BLOCK;		

			# save block ids
			set current_block_id_1 $b1
			set current_block_id_2 $b2
		
		
		} else {
		
			print_ERROR BEGIN_INTER_BLOCK;
			exit;
		
		}

	}

	proc END_INTER_BLOCK { } {
		variable level;
		variable level_name;

		if {[check_level 2 INTER_BLOCK] == 0 } {

			decr_level;			
						


		} else {
			
			print_ERROR END_INTER_BLOCK;
			exit;
		
		}

	}

# Master Equation File

	proc BEGIN_MEF { args } {
		variable level;
		variable level_name;
		variable lattice_cell;			
		variable lattice_vec_r1;
		variable lattice_vec_r2;
		variable lattice_vec_r3;
		
		
		if {[check_level 0 GROUND] == 0 } {
		
			incr_level MEF;		
			
			set lattice_cell(1) 0;
			set lattice_cell(2) 0;
			set lattice_cell(3) 0;
			set lattice_vec_r1(1) 0;
			set lattice_vec_r1(2) 0;
			set lattice_vec_r1(3) 0;
			set lattice_vec_r2(1) 0;
			set lattice_vec_r2(2) 0;
			set lattice_vec_r2(3) 0;
			set lattice_vec_r3(1) 0;
			set lattice_vec_r3(2) 0;
			set lattice_vec_r3(3) 0;
		
		
		} else {
		
			print_ERROR BEGIN_MEF;
			exit;
		
		}

	}

	proc END_MEF { } {
		variable level;
		variable level_name;

		if {[check_level 1 MEF] == 0 } {

			decr_level;			
						
		} else {
			
			print_ERROR END_MEF;
			exit;
		
		}
	}


	proc BEGIN_DM_INITIAL_CONDITION { b } {
		variable level;
		variable level_name;
		variable current_ic_id;
		variable dm_nr_states
		variable dm_elements;
		variable ic_count;

		if {[check_level 1 MEF] == 0 } {
		
			incr_level DM_INITIAL_CONDITION;		

			# save ic id
			set current_ic_id $b
			set dm_nr_states($current_ic_id) 0;
			incr ic_count;
		
		
		} else {
		
			print_ERROR BEGIN_DM_INITIAL_CONDITION;
			exit;
		
		}

	}

	proc END_DM_INITIAL_CONDITION { } {
		variable level;
		variable level_name;

		if {[check_level 2 DM_INITIAL_CONDITION] == 0 } {

			decr_level;			
						

		} else {
			
			print_ERROR END_DM_INITIAL_CONDITION;
			exit;
		
		}

	}

	proc BASIS_TYPE { type } {
		variable level;
		variable level_name;
		variable basis_type;
		variable current_ic_id;

		if {[check_level 2 DM_INITIAL_CONDITION] == 0 } {
		
			if {$type == "EXCITON"} {
		
				set basis_type($current_ic_id) 2
				
			} else if {$type == "EXCITON"} {
			
				set basis_type($current_ic_id) 3
			
			} else {
			
				print_ERROR "Unknown Basis type: $type";
				exit
			
			}

		} else {
			
			print_ERROR BASIS_TYPE;
			exit;
		
		}
	}


	proc NR_STATES { nr } {
		variable level;
		variable level_name;
		variable dm_nr_states;
		variable current_ic_id;
		variable dm_values_r;
		variable dm_values_i;

		if {[check_level 2 DM_INITIAL_CONDITION] == 0 } {
		
			set dm_nr_states($current_ic_id) $nr

			# emptying the dm_values
			for { set i 1 } { $i <= $nr } { incr i } {  
				for { set j 1 } { $j <= $nr } { incr j } {  			
					set dm_values_r($current_ic_id,$i,$j) 0;
					set dm_values_i($current_ic_id,$i,$j) 0;
				}
			}

		} else {
			
			print_ERROR NR_STATES;
			exit;
		
		}
	}


	proc ELEMENT { st1 st2 valr vali } {
		variable level;
		variable level_name;
		variable dm_values_r;
		variable dm_values_i;
		variable dm_nr_states;
		variable current_ic_id;

		if {[check_level 2 DM_INITIAL_CONDITION] == 0 } {
		
			if { $st1 > $dm_nr_states($current_ic_id) } {
			
				print_ERROR "Dimension of the DM in ELEMENT exceeded";
				exit;
		
			}

			if { $st2 > $dm_nr_states($current_ic_id) } {
			
				print_ERROR "Dimension of the DM in ELEMENT exceeded";
				exit;
		
			}

		
			set dm_values_r($current_ic_id,$st1,$st2) $valr
			set dm_values_i($current_ic_id,$st1,$st2) $vali

		} else {
			
			print_ERROR ELEMENT;
			exit;
		
		}
	}

	proc POP { st valr } {
		variable level_name;
		variable dm_values_r;
		variable dm_values_i;
		variable dm_nr_states;
		variable current_ic_id;

		if {[check_level 2 DM_INITIAL_CONDITION] == 0 } {

			ELEMENT $st $st $valr 0.0		

		} else {
			
			print_ERROR POP;
			exit;
		
		}
	
	}

	proc COH { st1 st2 valr vali } {
		variable level_name;
		variable dm_values_r;
		variable dm_values_i;
		variable dm_nr_states;
		variable current_ic_id;

		if {[check_level 2 DM_INITIAL_CONDITION] == 0 } {

			ELEMENT $st1 $st2 $valr $vali		

		} else {
			
			print_ERROR COH;
			exit;
		
		}
	
	}



#
#  NIS output
#

	proc putInteger { val } {

    	set ll [list s i $val];
		foreach a $ll {
			NOSE_conf::print $a;
		}	

	}

	proc putVector { en s1 } {
		upvar $en mat

		for {set i 1} {$i <= $s1 } {incr i} {
			append mm $mat($i) " ";
		}

		set ll [list v $s1 r $mm];
		foreach a $ll {
			NOSE_conf::print $a;
		}

	}

	proc putMatrix { en s1 s2 } {
		upvar $en mat

		for {set j 1} {$j <= $s2 } {incr j} {
			for {set i 1} {$i <= $s1 } {incr i} {
				append mm $mat($i,$j) " ";
			}
		}

		set ll [list m $s1 $s2 r $mm];
		foreach a $ll {
			NOSE_conf::print $a;
		}

	}

	proc putString { val } {

   		set l [string length $val]
    	set ll [list v $l c $val];
		foreach a $ll {
			NOSE_conf::print $a;
		}	

	}

	proc putEnergies { block_id } {	
		variable transition_energies;
		variable transition_ids;
		variable transition_count;
		
		set elist [list];
		
		for { set i 1 } { $i <= $transition_count($block_id) } { incr i } {
		
			set k [lindex $transition_ids($block_id) [expr $i - 1]];
			lappend elist $transition_energies($block_id,$k);
		
		}
		
		set ll [list v $transition_count($block_id) r $elist]
		foreach a $ll {
	    	NOSE_conf::print $a
	    	#puts $a;
		}
				
	}
	
	proc putGofts { block_id } {	
		variable transition_gofts;
		variable transition_ids;
		variable transition_count;
		
		set elist [list];
		
		for { set i 1 } { $i <= $transition_count($block_id) } { incr i } {
		
			set k [lindex $transition_ids($block_id) [expr $i - 1]];
			lappend elist $transition_gofts($block_id,$k);
		
		}
		
		set ll [list v $transition_count($block_id) i $elist]
		foreach a $ll {
	    	NOSE_conf::print $a
	    	#puts $a;
		}
				
	}
	
	proc putDisorders { block_id } {
	
		variable transition_disorders;
		variable transition_ids;
		variable transition_count;
		
		set elist [list];
		
		for { set i 1 } { $i <= $transition_count($block_id) } { incr i } {
		
			set k [lindex $transition_ids($block_id) [expr $i - 1]];
			lappend elist $transition_disorders($block_id,$k,$k);
		
		}
		
		set ll [list v $transition_count($block_id) i $elist]
		foreach a $ll {
	    	NOSE_conf::print $a
	    	#puts $a;
		}	
	
	}
	
	proc putRStrengths { block_id } {
		
		variable rstrength;
		variable transition_ids;
		variable transition_count;
		variable transition_ids;

		set j 0;
		foreach id $transition_ids($block_id) {
			incr j;
			
			if { [info exists rstrength($block_id,$id)] } { 
			
				set dd($j) $rstrength($block_id,$id);
				
			} else {
			
				set dd($j) 0.0;
				
			}
		
		    #puts $dd($j);
			
		}	
		
		
		putVector dd $j;
	
	}

	proc putDisWidths { block_id } {
		variable transition_ids;
		variable dis_width;
		variable transition_disorders;
	
			set j 0;
			foreach id1 $transition_ids($block_id) {
			    incr j;
			    set i 0;
				foreach id2 $transition_ids($block_id) {
					incr i;
					if { [info exists transition_disorders($block_id,$id1,$id2)] } {
	    				set dd($id1,$id2) $dis_width($transition_disorders($block_id,$id1,$id2));
	    			} else {
	    				set dd($id1,$id2) 0;
	    			}
				}		
			}	    
	    		
	    	putMatrix dd $j $i;
		
	
	}


	proc putBlockStart { } {

		set ll [list v 14 c SECTION_BLOCKS] 
		foreach a $ll {
			NOSE_conf::print $a;
		}

	}

	proc putBlockContinue { } {
	
		set ll [list v 14 c BLOCK_CONTINUE]
		foreach a $ll {
			NOSE_conf::print $a;
		}

	}

	proc putBlockEnd { } {

		set ll [list v 9 c BLOCK_END] 
		foreach a $ll {
			NOSE_conf::print $a;
		}

	}

	proc putPeriodicStart { } {

		set ll [list v 16 c SECTION_PERIODIC] 
		foreach a $ll {
			NOSE_conf::print $a;
		}

	}

	proc putPeriodicEnd { } {

		set ll [list v 12 c PERIODIC_END] 
		foreach a $ll {
			NOSE_conf::print $a;
		}

	}


	proc putInterBlockStart { } {

		set ll [list v 20 c SECTION_INTER_BLOCKS] 
		foreach a $ll {
			NOSE_conf::print $a;
			#puts $a;
		}

	}

	proc putInterBlockContinue { } {
	
		set ll [list v 20 c INTER_BLOCK_CONTINUE]
		foreach a $ll {
			NOSE_conf::print $a;
			#puts $a;
		}

	}

	proc putInterBlockEnd { } {

		set ll [list v 15 c INTER_BLOCK_END] 
		foreach a $ll {
			NOSE_conf::print $a;
			#puts $a;
		}

	}



	#
	# Writes the block id into NIS
	#
	proc putBlockId { i } {

		putInteger $i

	}
	
	#
	# Writes the block id into NIS
	#
	proc putInterBlockId { i j } {

		putInteger $i
		putInteger $j
		#puts $i
		#puts $j

	}

	#
	# Writes the goft id into NIS
	#
	proc putGoftId { i } {

		putInteger $i

	}

	#
	# Writes the value determining the presence of the PERIODIC section
	#
	proc putPeriodic { i } {
	
		putInteger $i;
	
	}

	#
	# Writes the dimesion of the peridic structure
	#
	proc putDimension { i } {
	
		putInteger $i;
	
	}

	#
	# Writes the dimesion of the peridic structure
	#
	proc putPeriodicity { i } {
	
		putInteger $i;
	
	}



	#
	#
	#
	proc putCoords { block_id } {
		variable transition_coordinates;
		variable transition_ids;
		variable transition_count;	

		for { set i 1 } { $i <= 3} { incr i } {
			
			set j 0;
			foreach id $transition_ids($block_id) {

			    incr j;
	    		set dd($j,$i) $transition_coordinates($block_id,$id,$i);
						
			}	    
	    		
	    	
	    }    

	    putMatrix dd $j 3;
	
	}
	
	
	proc putHarmonicOscillators { block_id } {
	
		variable transition_vibrational_levels;
		variable transition_ids;
		variable transition_count;
		variable qo_levels;
		variable qo_frequency;
		variable qo_reorganization_energy;
		variable qo_huang_rhys_factor;
		
		set elist [list];
		
		set qo_levels(0) 0;
		set qo_frequency(0) 0;
		set qo_reorganization_energy(0) 0;
		set qo_huang_rhys_factor(0) 0;
		
		for { set i 1 } { $i <= $transition_count($block_id) } { incr i } {
		
			set k [lindex $transition_ids($block_id) [expr $i - 1]];
			set kk $transition_vibrational_levels($block_id,$k);

			if { [info exists qo_reorganization_energy($kk) ] && $qo_frequency($kk) != 0  } {

				set qo_huang_rhys_factor($kk) [ expr $qo_reorganization_energy($kk) / $qo_frequency($kk) ];
				
			} elseif { [info exists qo_huang_rhys_factor($kk) ] && $qo_frequency($kk) != 0 }  { 

				set qo_reorganization_energy($kk) [ expr $qo_huang_rhys_factor($kk) * $qo_frequency($kk) ];
				
			} else {
				set qo_reorganization_energy($kk) 0
				set qo_huang_rhys_factor($kk) 0
			}

			lappend elist $qo_levels($kk) $qo_frequency($kk) $qo_reorganization_energy($kk) $qo_huang_rhys_factor($kk);
		
		}
		
		set ll [list v [expr $transition_count($block_id) * 4] r $elist]
		foreach a $ll {
		    NOSE_conf::print $a
		    #puts $a;
		}	
	
	}
		
	
	#
	# Lattice vectors
	#
	proc putRR { } {
	
		variable lattice_vec_r1;
		variable lattice_vec_r2;
		variable lattice_vec_r3;
		
		putVector lattice_vec_r1 3;
		putVector lattice_vec_r2 3;
		putVector lattice_vec_r3 3;
	
	}

	#
	# Neighbour cell coupling 1D
	#
	proc putNCCoupling1D { } {
	
		variable cell_coupling;
		variable transition_count;

		
		for { set i 1} {$i <= $transition_count(1)} { incr i} {
		
			for { set j 1 } { $j <= $transition_count(1)} {incr j} {
			
				set mat($i,$j) $cell_coupling(1,0,$i,$j);
				puts "$cell_coupling(1,0,$i,$j);";
			
			}
			
		}
		

		putMatrix mat $transition_count(1) $transition_count(1);
		
	
	}

	#
	# Neighbour cell coupling 1D
	#
	proc putNCCoupling2D { } {
	
		variable cell_coupling;
		variable transition_count;

		
		for { set i 1} {$i <= $transition_count(1)} { incr i} {
		
			for { set j 1 } { $j <= $transition_count(1)} {incr j} {
			
				set mat($i,$j) $cell_coupling(1,0,$i,$j);
				puts "$cell_coupling(1,0,$i,$j);";
			
			}
			
		}
		

		putMatrix mat $transition_count(1) $transition_count(1);
		
		for { set i 1} {$i <= $transition_count(1)} { incr i} {
		
			for { set j 1 } { $j <= $transition_count(1)} {incr j} {
			
				set mat($i,$j) $cell_coupling(0,1,$i,$j);
				puts "$cell_coupling(0,1,$i,$j);";
			
			}
			
		}
		

		putMatrix mat $transition_count(1) $transition_count(1);

		for { set i 1} {$i <= $transition_count(1)} { incr i} {
		
			for { set j 1 } { $j <= $transition_count(1)} {incr j} {
			
				set mat($i,$j) $cell_coupling(1,1,$i,$j);
				puts "$cell_coupling(1,1,$i,$j);";
			
			}
			
		}
		

		putMatrix mat $transition_count(1) $transition_count(1);
	

		for { set i 1} {$i <= $transition_count(1)} { incr i} {
		
			for { set j 1 } { $j <= $transition_count(1)} {incr j} {
			
				set mat($i,$j) $cell_coupling(-1,1,$i,$j);
				puts "$cell_coupling(-1,1,$i,$j);";
			
			}
			
		}
		

		putMatrix mat $transition_count(1) $transition_count(1);
	
	}


	#
	#
	#
	proc putDipole { block_id } {
		variable transition_dipoles;
		variable transition_ids;
		variable transition_count;	

		
		for { set i 1 } { $i <= 3} { incr i } {
			
			set j 0;
			foreach id $transition_ids($block_id) {

			    incr j;
	    		set dd($j) $transition_dipoles($block_id,$id,$i);
						
			}	    
	    		
	    	putVector dd $j;
	    	
	    }    
		
	
	}

	proc putDLengths { block_id } {	
		variable transition_dlengths;
		variable transition_ids;
		variable transition_count;
		
		set elist [list];
		
		for { set i 1 } { $i <= $transition_count($block_id) } { incr i } {
		
			set k [lindex $transition_ids($block_id) [expr $i - 1]];
			lappend elist $transition_dlengths($block_id,$k);
		
		}
		
		set ll [list v $transition_count($block_id) r $elist]
		foreach a $ll {
	    	NOSE_conf::print $a
	    	# puts $a;
		}
				
	}
	

	#
	# Writes coupling matrix into NIS 
	#
	proc putCouplings { block_id } {
		variable transition_ids;
		variable transition_count;	
		variable couplings;

		for { set a 1 } { $a <= $transition_count($block_id) } { incr a } {
			for { set b 1 } { $b <= $transition_count($block_id) } { incr b } {
		
				if { [info exists couplings($block_id,$a,$b)] } {
					set mat($a,$b) $couplings($block_id,$a,$b);
				} else {
					set mat($a,$b) 0.0;
				}
			
			}
			
		}

		putMatrix mat $transition_count($block_id) $transition_count($block_id); 

	}

	#
	# Writes inter coupling matrix into NIS 
	#
	proc putInterCouplings { block_id1 block_id2 } {
		variable transition_ids;
		variable transition_count;	
		variable int_couplings;

		for { set a 1 } { $a <= $transition_count($block_id1) } { incr a } {
			for { set b 1 } { $b <= $transition_count($block_id2) } { incr b } { 
		
				#puts $block_id1
				#puts $block_id2
				
				if { [info exists int_couplings($block_id1,$block_id2,$a,$b)] } {
					set mat($a,$b) $int_couplings($block_id1,$block_id2,$a,$b);
					#puts "Exists"
				} else {
					set mat($a,$b) 0.0;
				}
			
			}
			
		}

		putMatrix mat $transition_count($block_id1) $transition_count($block_id2); 
		#puts "$transition_count($block_id1) $transition_count($block_id2) ";
	}
	

	proc putGoftStart { } {

		putString SECTION_GOFTS; 

	}

	proc putGoftContinue { } {
	
		putString GOFT_CONTINUE;

	}

	proc putGoftEnd { } {

		putString GOFT_END; 

	}

	proc putInitCondStart { } {

		putString SECTION_DM_INITIAL_CONDITION; 

	}

	proc putInitCondEnd { } {

		putString DM_INITIAL_CONDITION_END; 

	}

	proc putBasisType { i } {
	
		putInteger $i;
	
	}

	proc putNrStates { i } {
	
		putInteger $i;
	
	}
	
	
	proc putInitialCondition { ic_id } {	
		variable dm_nr_states;
		variable dm_values_r;
		variable dm_values_i;
		
		set elist [list];
		
		for { set i 1 } { $i <= $dm_nr_states($ic_id) } { incr i } {
		
			lappend elist $dm_values_r($ic_id,$i,$i);
			lappend elist2 $dm_values_i($ic_id,$i,$i);
		
		}
		
		set ll [list v $dm_nr_states($ic_id) r $elist]
		set ll2 [list v $dm_nr_states($ic_id) r $elist2]
		foreach a $ll {
	    	NOSE_conf::print $a
	    	#puts $a;
		}
		foreach a $ll2 {
	    	NOSE_conf::print $a
	    	#puts $a;
		}
				
	}	


	proc putGoftMode { args } {
	
		set state type;
		set tp    "";
		set count 0;
		
		foreach arg $args {
		
			switch -- $state {
			
				type {  
					
					set tp $arg
					set state param
				
				}
				
				param {  
				
					switch -- $tp {
					
						BROWNIAN {
							set lpar $arg;
							if { [llength $lpar] == 2 } {
								set par(1) [lindex $lpar 0];
								set par(2) [lindex $lpar 1];									
								putBrownianMode $par(1) $par(2);
								set state type;
							}	
						}

						BROWNIAN_NO_MATSUBARA {
							set lpar $arg;
							if { [llength $lpar] == 2 } {
								set par(1) [lindex $lpar 0];
								set par(2) [lindex $lpar 1];									
								putBrownianNoMatsubaraMode $par(1) $par(2);
								set state type;
							}	
						}

						BROWNIAN_LOW_TEMP_HIERARCHY {
							set lpar $arg;
							if { [llength $lpar] == 2 } {
								set par(1) [lindex $lpar 0];
								set par(2) [lindex $lpar 1];									
								putBrownianLowTempHierarchyMode $par(1) $par(2);
								set state type;
							}	
						}
					
						GAUSSIAN {
							set lpar $arg;
							if { [llength $lpar] == 2 } {
								set par(1) [lindex $lpar 0];
								set par(2) [lindex $lpar 1];									
								putGaussianMode $par(1) $par(2);
								set state type;
							}	
						}
						
						BROWNIAN_GENERAL {
							set lpar $arg;
							if { [llength $lpar] == 3 } {
								set par(1) [lindex $lpar 0];
								set par(2) [lindex $lpar 1];
								set par(3) [lindex $lpar 2];
								putGeneralBrownianMode $par(1) $par(2) $par(3);
								set state type;
							}	
						}
						
						BROWNIAN_UNDERDAMPED {
							set lpar $arg;
							if { [llength $lpar] == 2 } {
								set par(1) [lindex $lpar 0];
								set par(2) [lindex $lpar 1];
								putUnderdampedBrownianMode $par(1) $par(2);
								set state type;
							}	
						}						

						DELTA {
							set lpar $arg;
							if { [llength $lpar] == 1 } {
								set par(1) [lindex $lpar 0];
								putDeltaMode $par(1);
								set state type;
							}	

						}
					
					}
				
				}
			
			}
		
		}
	
	}

	proc putModeEnd { } {
	
		putString MODE_END;
	
	}


	proc putEndSections { } {

		putString END_SECTIONS;

	} 


	proc putBrownianMode { reorg tc } {

		putString BROWNIAN;
	
		set ll [list v 2 r [list $reorg $tc]]
		foreach a $ll {
			NOSE_conf::print $a;
		}	

	}
	
	proc putBrownianNoMatsubaraMode { reorg tc } {

		putString BROWNIAN_NO_MATSUBARA;
	
		set ll [list v 2 r [list $reorg $tc]]
		foreach a $ll {
			NOSE_conf::print $a;
		}	

	}

	proc putBrownianLowTempHierarchyMode { reorg tc } {

		putString BROWNIAN_LOW_TEMP_HIERARCHY;
	
		set ll [list v 2 r [list $reorg $tc]]
		foreach a $ll {
			NOSE_conf::print $a;
		}	

	}

	proc putGaussianMode { reorg tc } {
	
		putString GAUSSIAN;
	
		set ll [list v 2 r [list $reorg $tc]]
		foreach a $ll {
			NOSE_conf::print $a;
		}		
	
	}
	
	proc putGeneralBrownianMode { reorg tc omeg } {

		putString BROWNIAN_GENERAL;
	
		set ll [list v 3 r [list $reorg $tc $omeg]]
		foreach a $ll {
			NOSE_conf::print $a;
		}	

	}	
	
	proc putUnderdampedBrownianMode { reorg omeg } {

		putString BROWNIAN_UNDERDAMPED;
	
		set ll [list v 2 r [list $reorg $omeg]]
		foreach a $ll {
			NOSE_conf::print $a;
		}	

	}

	proc putDeltaMode { gm } {

		putString DELTA;
	
		set ll [list v 1 r [list $gm]]
		foreach a $ll {
			NOSE_conf::print $a;
		}	

	}




	#
	# Dipole - dipole interaction energy
	#
	proc dipole_dipole { tr1 tr2 } {
		variable transition_coordinates;
		variable transition_dipoles;
		variable transition_dlengths;
		variable current_block_id;
		variable nose_home;
		variable nexdir;
		variable lattice_cell;
		variable periodicity_dim;
		variable lattice_vec_r1;
		variable lattice_vec_r2;
		variable lattice_vec_r3;
		
		global   use_expect_package;
		
		
		for { set i 1 } { $i <= 3 } { incr i } {

			set r1($i) $transition_coordinates($current_block_id,$tr1,$i);
			set r2($i) $transition_coordinates($current_block_id,$tr2,$i);
			
			
			set x1 [expr ( $lattice_vec_r1($i) * $lattice_cell(1) )];
			if {$periodicity_dim > 1} {
				set x1 [expr $x1 + ( $lattice_vec_r2($i) * $lattice_cell(2) ) ];
				if {$periodicity_dim > 2} {
					set x1 [expr $x1 + ( $lattice_vec_r3($i) * $lattice_cell(3) ) ];
				}	
			}	
#			puts "x1 = $x1";
				
			set r2($i) [expr $r2($i) + $x1];
				

						
		}
		
		for { set i 1 } { $i <= 3 } { incr i } {
			set d1($i) $transition_dipoles($current_block_id,$tr1,$i);
			set d2($i) $transition_dipoles($current_block_id,$tr2,$i);			
		}
		
		set l1 $transition_dlengths($current_block_id,$tr1);
		set l2 $transition_dlengths($current_block_id,$tr2);
		
		set n1 [expr sqrt( ( $d1(1) * $d1(1) ) + ( $d1(2) * $d1(2) ) + ( $d1(3) * $d1(3) ) ) ];
		set n2 [expr sqrt( ( $d2(1) * $d2(1) ) + ( $d2(2) * $d2(2) ) + ( $d2(3) * $d2(3) ) ) ];
		 
		set f1 [expr $l1 / $n1];
		set f2 [expr $l2 / $n2];
		
		
		for { set i 1 } { $i <= 3 } { incr i } {
			set ddip1($i) [expr $d1($i) * $f1 ];
			set ddip2($i) [expr $d2($i) * $f2 ];
			set dR($i) [expr $r1($i) - $r2($i) ];
		}
		
		set RRnorm [expr sqrt( ( $dR(1) * $dR(1) ) + ( $dR(2) * $dR(2) ) + ( $dR(3) * $dR(3) ) ) ];
		
		for { set i 1 } { $i <= 3 } { incr i } {
			set dR($i) [expr $dR($i) / $RRnorm ];
		}
		
		set dprodd1d2 [expr  ( $ddip1(1) * $ddip2(1) ) + ( $ddip1(2) * $ddip2(2) ) + ( $ddip1(3) * $ddip2(3) ) ];
		set dprodd1r0 [expr  ( $ddip1(1) * $dR(1) ) + ( $ddip1(2) * $dR(2) ) + ( $ddip1(3) * $dR(3) ) ];
		set dprodd2r0 [expr  ( $ddip2(1) * $dR(1) ) + ( $ddip2(2) * $dR(2) ) + ( $ddip2(3) * $dR(3) ) ];
		
		set ret1 0;
		set ret2 [expr  5034.1153 * ($dprodd1d2 - 3 * $dprodd1r0 * $dprodd2r0) / ($RRnorm * $RRnorm * $RRnorm) ];
		
		set ex [file join $nexdir dip.x];
		
#		# There is problem sometimes if they are defined only in local scopes of if under this sign
#		# so I define them here somehow. They are overwritten in any case.
#		set ret1 42
#		set ret2 42
#				
#		if { $use_expect_package == 1 } {		
#			if {[catch {
#			set timeout 1
#			spawn -noecho $ex;
#			log_user 0
#		
#			for { set i 1 } { $i <= 3 } { incr i } {
#				send -- $r1($i);
#				send "\n";
#				send -- $r2($i);
#				send "\n";
#				send -- [expr $d1($i) * $f1 ];
#				send "\n";
#				send -- [expr $d2($i) * $f2 ];
#				send "\n";
#			}
#			expect -re "err\[^\n]*(-?\[0-9]\[0-9]*\.?\[0-9]*e?-?\[0-9]*)";
#			set ret1 [string replace $expect_out(0,string) 0 4];
#			
#			expect -re "vvv\[^\n]*(-?\[0-9]\[0-9]*\.?\[0-9]*e?-?\[0-9]*)";
#			set ret2 [string replace $expect_out(0,string) 0 4]
#			
#			expect
#			log_user 1
#			} result ] } {
#			puts "error in expect section: $result"
#			}
#		} else {
#			set res [open "| $ex " w+];
#			puts $ex;
#			for { set i 1 } { $i <= 3 } { incr i } {
#				puts $res $r1($i);
#				puts $res $r2($i);			
#				puts $res [expr $d1($i) * $f1 ];
#				puts $res [expr $d2($i) * $f2 ];				
#			}
#			flush $res;
#			
#			puts "THERE IS A TROUBLE HERE!";
#			puts "Edit nose/src/scripts/nose.tcl and see the comment or do not use dipole calculation in your config-file in no-Expect mode.";
#			exit;
#			
#			# This part doesn't work now. ret1 and ret2 should			
#			# contain only number, but dip.x routine returns 
#			# also some characters which help Expect to identify 
#			# returned numbers. Only thing that needs to be done
#			# is to extract the numbers from ret1 and ret2. 
#			# I didn't do it so far, because pipes on ginny are
#			# broken and I cannot properly test it.
#			 
#			set ret1 [gets $res];
#			set ret2 [gets $res];
#					
#			close $res;	
#			puts "DONE!";
#		}
		
		if { $ret1 != 0 } {
		
			puts "Error in dipole_dipole calculation";
		
		} 
			    
		return $ret2;

	}



	namespace export BEGIN_SSF;
	namespace export END_SSF;
	namespace export BEGIN_BLOCK;
	namespace export END_BLOCK;
	namespace export TRANSITION;
	namespace export ROTATIONAL_STRENGTH;
	namespace export COUPLING;
	namespace export UNITS;
	namespace export UNITS_ENERGY;
	namespace export UNITS_TIME;
	namespace export BEGIN_CORRF;
	namespace export END_CORRF;
	namespace export BEGIN_MODE;
	namespace export END_MODE;
	namespace export MODE_TYPE;
	namespace export REORGANIZATION_ENERGY;
	namespace export CORRELATION_TIME;
	namespace export BROWNIAN_OMEGA;
	namespace export GAMMA;
	namespace export BEGIN_DISORDER;
	namespace export BEGIN_QUANTUM_OSCILLATOR
	namespace export END_QUANTUM_OSCILLATOR
	namespace export QO_HUANG_RHYS_FACTOR
	namespace export QO_OSCILLATOR
	namespace export QO_REORGANIZATION_ENERGY
	namespace export QO_FREQUENCY
	namespace export QO_LEVELS	
	namespace export END_DISORDER; 
	namespace export TYPE;
	namespace export TARGETS;
	namespace export WIDTH;
	namespace export BEGIN_GROUP;
	namespace export END_GROUP;
	namespace export MEMBERS;
	namespace export DISORDER;
	namespace export BEGIN_INTER_BLOCK;
	namespace export END_INTER_BLOCK;
	namespace export BEGIN_PERIODIC;
	namespace export END_PERIODIC;
	namespace export DIMENSIONS;
	namespace export PERIODICITY;
	namespace export R1;
	namespace export R2;
	namespace export R3;
	namespace export CELL;
	namespace export TURN;
	
	namespace export BEGIN_MEF;
	namespace export END_MEF;
	namespace export BEGIN_DM_INITIAL_CONDITION;
	namespace export END_DM_INITIAL_CONDITION;
	namespace export BASIS_TYPE;
	namespace export NR_STATES;
	namespace export ELEMENT;
	namespace export POP;
	namespace export COH;
	
	namespace export print_LEVEL;
	 
	namespace export dipole_dipole;

} ;# END namespace NOSE_ssf


global env

##############################################
#  IMPORT NAMES
##############################################

namespace import NOSE_conf::moduleName;
namespace import NOSE_conf::moduleInpFile;
namespace import NOSE_conf::inputFile;
namespace import NOSE_conf::gridSteps;
namespace import NOSE_conf::gridExtent;
namespace import NOSE_conf::extendedGridExtent;
namespace import NOSE_conf::parallel;
namespace import NOSE_conf::tau;
namespace import NOSE_conf::dgamma;
namespace import NOSE_conf::temperature;
namespace import NOSE_conf::timeStep;
namespace import NOSE_conf::logLevel;
namespace import NOSE_conf::pulsePolarization
namespace import NOSE_conf::detectionPolarization;
namespace import NOSE_conf::outputName;
namespace import NOSE_conf::outputDir;
namespace import NOSE_conf::output;
namespace import NOSE_conf::outLog;
namespace import NOSE_conf::outLogIds;
namespace import NOSE_conf::saveNIS;
namespace import NOSE_conf::saveGoft;
namespace import NOSE_conf::restartFreq;
namespace import NOSE_conf::realizations;
namespace import NOSE_conf::tempDir;
namespace import NOSE_conf::deleteTemp;
namespace import NOSE_conf::units;
namespace import NOSE_conf::pulseFrequency;
namespace import NOSE_conf::putParametricRecordName;
namespace import NOSE_conf::rwa;
namespace import NOSE_conf::moduleMethod;
namespace import NOSE_conf::completeResultFile;
namespace import NOSE_conf::parallel_options;
namespace import NOSE_conf::excPosition;
namespace import NOSE_conf::molecule;
namespace import NOSE_conf::pathway;
namespace import NOSE_conf::runs;
namespace import NOSE_conf::spec_wini-dw-wst;
namespace import NOSE_conf::localBasis;
namespace import NOSE_conf::relaxation;
namespace import NOSE_conf::secular;
namespace import NOSE_conf::dephasing;
namespace import NOSE_conf::feeding;
namespace import NOSE_conf::draining;

namespace import NOSE_ssf::BEGIN_SSF;
namespace import NOSE_ssf::END_SSF;
namespace import NOSE_ssf::BEGIN_BLOCK;
namespace import NOSE_ssf::END_BLOCK;
namespace import NOSE_ssf::TRANSITION;
namespace import NOSE_ssf::ROTATIONAL_STRENGTH;
namespace import NOSE_ssf::COUPLING;
namespace import NOSE_ssf::UNITS;
namespace import NOSE_ssf::UNITS_ENERGY;
namespace import NOSE_ssf::UNITS_TIME;
namespace import NOSE_ssf::BEGIN_CORRF;
namespace import NOSE_ssf::END_CORRF;
namespace import NOSE_ssf::MODE_TYPE;
namespace import NOSE_ssf::CORRELATION_TIME;
namespace import NOSE_ssf::REORGANIZATION_ENERGY;
namespace import NOSE_ssf::BROWNIAN_OMEGA;
namespace import NOSE_ssf::BEGIN_DISORDER;
namespace import NOSE_ssf::END_DISORDER;
namespace import NOSE_ssf::BEGIN_QUANTUM_OSCILLATOR
namespace import NOSE_ssf::END_QUANTUM_OSCILLATOR
namespace import NOSE_ssf::QO_HUANG_RHYS_FACTOR
namespace import NOSE_ssf::QO_OSCILLATOR
namespace import NOSE_ssf::QO_REORGANIZATION_ENERGY
namespace import NOSE_ssf::QO_FREQUENCY
namespace import NOSE_ssf::QO_LEVELS
namespace import NOSE_ssf::TYPE;
namespace import NOSE_ssf::TARGETS;
namespace import NOSE_ssf::WIDTH;
namespace import NOSE_ssf::BEGIN_GROUP;
namespace import NOSE_ssf::END_GROUP;
namespace import NOSE_ssf::MEMBERS;
namespace import NOSE_ssf::DISORDER;
namespace import NOSE_ssf::BEGIN_INTER_BLOCK;
namespace import NOSE_ssf::END_INTER_BLOCK;
namespace import NOSE_ssf::BEGIN_MODE;
namespace import NOSE_ssf::END_MODE;
namespace import NOSE_ssf::GAMMA;

namespace import NOSE_ssf::BEGIN_PERIODIC;
namespace import NOSE_ssf::END_PERIODIC;
namespace import NOSE_ssf::DIMENSIONS;
namespace import NOSE_ssf::PERIODICITY;
namespace import NOSE_ssf::R1;
namespace import NOSE_ssf::R2;
namespace import NOSE_ssf::R3;
namespace import NOSE_ssf::CELL;
namespace import NOSE_ssf::TURN;

namespace import NOSE_ssf::BEGIN_MEF;
namespace import NOSE_ssf::END_MEF;
namespace import NOSE_ssf::BEGIN_DM_INITIAL_CONDITION;
namespace import NOSE_ssf::END_DM_INITIAL_CONDITION;
namespace import NOSE_ssf::BASIS_TYPE;
namespace import NOSE_ssf::NR_STATES;
namespace import NOSE_ssf::ELEMENT;
namespace import NOSE_ssf::POP;
namespace import NOSE_ssf::COH;

namespace import NOSE_ssf::print_LEVEL;


namespace import NOSE_ssf::dipole_dipole;

global nos 

set nos 0;
 
proc Reader { pipe } {
	global done
	global nos
	if {[eof $pipe]} {
		if {[catch {close $pipe} result]} {
			puts stderr $result
		}
		set done 1
		return
	}
	gets $pipe line;
	if { [string compare $line "<!--External"] == 0 } {
		set nos 1	
	} 
	if { [string compare $line "External--!>"] == 0 } {
		set nos 2
	}
	if { $nos == 0 } {
		puts $line;
	} else {
		puts " b: $line";
		if { $nos == 2 } {
			set nos 0;
		}
	}
}


##############################################
#  Default configuration
##############################################

output default_log
#completeResultFile no

#set genmat 0; 


##############################################
#  Use the configuration file
##############################################


#
# Put out help and exit
#
if { [llength $argv] == 0 } {

   puts " ";
   puts "Usage: nose \[options\] \[configfile\]";
   puts "Report bugs to: tomas.mancal@mff.cuni.cz";
   puts " ";
   puts "Type: nose --help   to get more information"
   puts " ";  
   exit  

}

set confile [lindex $argv 0];

if { [string compare $confile "--help"] == 0 } {

   puts " "
   puts "Usage: nose \[options\] \[configfile\]";
   puts "Report bugs to: tomas.mancal@mff.cuni.cz";
   puts " ";
   puts "Options:";
   puts " --help ................... prints this message";
   puts " --info ................... prints a more detailed message";
   puts " --tutorial \[dirname\] ..... generates a tutorial directory";
   puts "                            default dirname is \"tutorial\"";
   puts " --debugging configfile ... runs nose in debugging mode";
   puts " ";
   exit 

}

if { [string compare $confile "--info"] == 0 } {

   puts " "
   puts "Usage: nose \[options\] \[configfile\]";
   puts "Report bugs to: tomas.mancal@mff.cuni.cz";
   puts " ";
   puts "Options:";
   puts " --help ................... prints this message";
   puts " --info ................... prints a more detailed message";
   puts " --tutorial \[dirname\] ..... generates a tutorial directory";
   puts "                            default dirname is \"tutorial\"";
   puts " --debugging configfile ... runs nose in debugging mode";
   puts " ";
   puts "Info: ";
   puts " Compilation options : $fcflags";
   puts " Using mpi           : $mpi";
   puts " ";
   exit 

}



#
# Generate a tutorial directory and exit
#
if { [string compare $confile "--tutorial"] == 0 } {

	if { [llength $argv] == 1 } {
	
		set tdir [file join ./ tutorial] 
			
	}   
	
	if { [llength $argv] >= 2 } {
	
		set tdir [lindex $argv 1];
		set tdir [file join ./ $tdir];
		
	} 

	if { [file exists $tdir] == 1 } {
	
		puts "Error: The intended tutorial directory already exists!"
		puts "Choose a different name or delete the existing directory."
		exit
	
	}
	puts "Generating tutorial directory $tdir ";

	set scrname [info script];
	set dname [file join [file dirname $scrname] ../tutorial];
	
	file copy $dname $tdir
	
	puts $dname
	
	exit;
	
}

#
# Sets the nose mode to debugging
#
set dbg no;
if { [string compare $confile "--debugging"] == 0 } {

	if { [llength $argv] == 1 } {
	
		puts " ";
		puts "Error: this option has to be followed by configuration file name";
		puts "Type: nose --help";
		puts " ";
		exit; 
					
	} 
	
	if { [llength $argv] >= 2 } {
	
		set confile [lindex $argv 1];

		#
		# Test for mpi
		#
		if { [string compare $mpi yes] == 0 } {		
			puts " ";
			puts "Warning: Debugging of parallel nose is not supported";
			puts "Reconfigure with \"--with-mpi=no\" option, and recompile ";
			puts " ";
			#exit;
		}
		
		#
		# Test for -g option
		#
		puts $fcflags
		if { [string match "*-g*" $fcflags] == 1 } {
		
		
		} else {
			puts " ";
			puts "Error: For debugging, NOSE has to be compiled with -g option only";
			puts "Reconfigure with \"--enable-debugging\" option, and recompile ";
			puts " ";
			exit;		
		}
		
	}

	set dbg yes;
	
} 

#
# Checking for testing mode
#
if { [string compare $confile "--check"] == 0 } {

    set confile [lindex $argv 1];

    set exdir [file join .. .. src pert];

}


#
# Just save nis to a file
#
set nis_only no
if { [string compare $confile "--nis-file"] == 0 } {

	set nis_only yes;
	
    set confile [lindex $argv 1];	

}


#
# Source the conf file
#

puts "Info   : Output from the NOSE configuration file follows ... "
puts "->------------------------------------------------------------------------------"
switch  [catch {
	source $confile;
} result ] {
	1 {	puts stderr $result;}
	default {               }
}
puts "-<------------------------------------------------------------------------------"
puts "Info   : ... configuration file finished"

##############################################
#  Source things from the conf file
##############################################
if { $NOSE_conf::smth2source == 1 } {
	source $NOSE_conf::to_source;
}


################################################
# Parse input file
################################################

# defaults in the conf file

set NOSE_conf::e_units wn

#################################################
# Input file pre-defaults
#################################################
#
# Source the input file
#
puts "Info   : Output from the NOSE input file follows ... "
puts "->------------------------------------------------------------------------------"
source $NOSE_conf::input_file;
puts "-<------------------------------------------------------------------------------"
puts "Info   : ... input file finished"

#################################################
# Input file post-defaults
#################################################

BEGIN_SSF

BEGIN_DISORDER 0
TYPE GAUSSIAN
WIDTH 0
END_DISORDER
 
END_SSF

##############################################
# executable and options
##############################################
set nose_home         $nhome ;# $env(NOSE_HOME); 

set nose_bin          $exdir; # [file join $nose_home bin];

if {$NOSE_conf::parallel_execution} {

	set mpi yes

} else {

	set mpi no

}

if { [string compare $dbg yes] == 0 } {

	set driver             ""
	set driver_options     ""
	set debug_exec [file join $nose_bin ndrv.x];

} else {

	if { [string compare $nis_only yes] == 0 } {
	
		set driver   ""
		set driver_options ""
	
	} else {

		set driver             ndrv.x
		set driver_options     ""
		
	}
	
}	



if { [string compare $mpi yes] == 0 } {
 
	set loader             mpirun 
    set loader_options     $NOSE_conf::conf_parallel_options

	
	
} else {

	set loader      ""
	set loader_options ""

}


#################################################
#
#################################################


##############################################
#  Start nose driver
##############################################
if {[string compare $driver ""] != 0} {

set driver [file join $nose_bin $driver]; 

}

if { [string compare $dbg no ] == 0} {

	puts "Info   : Debugging executable  $loader $loader_options $driver $driver_options";

}

######## log directory ##############

if { ![file exists log] } {  
   file mkdir log;
}

##############################################
# Stream output
##############################################

set out [open nis.in w ];
set nis_in nis.in;


#########################################################################################3
#  NIS start
###############################################################################

NOSE_conf::putInit out

####################################################################
# Configuration
####################################################################

NOSE_conf::putRecord module_name
NOSE_conf::putRecord temp
NOSE_conf::putRecord timeStep
NOSE_conf::putRecord grid_steps
NOSE_conf::putRecord grid_extent
NOSE_conf::putRecord saveGoft
NOSE_conf::putRecord saveNIS
NOSE_conf::putRecord parallel
NOSE_conf::putRecord inp_file
NOSE_conf::putRecord module_inp_file
NOSE_conf::putRecord logLevel
NOSE_conf::putRecord output_dir
#NOSE_conf::putRecord completeResultFile
NOSE_conf::putRecord realizations
NOSE_conf::putRecord restartFreq

NOSE_conf::putParametricRecord output
NOSE_conf::putRecord tau_projector

NOSE_conf::putRecord rwa

#########################################################################################
# SSF File
#########################################################################################

# blocks
NOSE_ssf::putBlockStart
for { set i 1 } { $i <= $NOSE_ssf::block_count } { incr i } {
	NOSE_ssf::putBlockContinue
	NOSE_ssf::putBlockId $i 
	NOSE_ssf::putEnergies $i
	NOSE_ssf::putCoords $i
	NOSE_ssf::putDipole $i
	NOSE_ssf::putDLengths $i                 ;# !!!!!!!!!!!!!!!!!!!!!!!!!1
	NOSE_ssf::putCouplings $i
	NOSE_ssf::putGofts $i
	NOSE_ssf::putDisorders $i
	NOSE_ssf::putDisWidths $i
	NOSE_ssf::putRStrengths $i
	NOSE_ssf::putHarmonicOscillators $i
}
NOSE_ssf::putBlockEnd

NOSE_ssf::putPeriodicStart
NOSE_ssf::putPeriodic $NOSE_ssf::periodic;
if { ($NOSE_ssf::block_count == 1) && ($NOSE_ssf::periodic == 1) } {
	NOSE_ssf::putDimension $NOSE_ssf::periodicity_dim;
	NOSE_ssf::putPeriodicity $NOSE_ssf::periodicity;
	NOSE_ssf::putRR;  
	
	if {$NOSE_ssf::periodicity_dim == 1} {
	
		NOSE_ssf::putNCCoupling1D;
		
	} elseif {$NOSE_ssf::periodicity_dim == 2} {
	
		NOSE_ssf::putNCCoupling2D;
	
	}

}
NOSE_ssf::putPeriodicEnd

# inter blocks
NOSE_ssf::putInterBlockStart
for { set i 1 } { $i <= $NOSE_ssf::block_count } { incr i } {
for { set j 1 } { $j <= $NOSE_ssf::block_count } { incr j } {
	if { $i != $j } {
		NOSE_ssf::putInterBlockContinue
		NOSE_ssf::putInterBlockId $i $j 
		NOSE_ssf::putInterCouplings $i $j
	}
}
}
NOSE_ssf::putInterBlockEnd

set NOSE_conf::e_units wn;

# gofts
NOSE_ssf::putGoftStart
for { set i 1 } { $i <= $NOSE_ssf::goft_count } { incr i } {
	NOSE_ssf::putGoftContinue
	NOSE_ssf::putGoftId $i
	for {set j 1 } { $j <= $NOSE_ssf::mcounts($i) } {incr j} {
	   
		switch -glob -- $NOSE_ssf::mode_type($i,$j) {

			BROWNIAN_GENERAL { set parms [list $NOSE_ssf::reorganization_energy($i,$j)  $NOSE_ssf::correlation_time($i,$j)  $NOSE_ssf::brownian_omega($i,$j)];
					   NOSE_ssf::putGoftMode BROWNIAN_GENERAL $parms; }

			BROWNIAN { set parms [list $NOSE_ssf::reorganization_energy($i,$j)  $NOSE_ssf::correlation_time($i,$j)];
					   NOSE_ssf::putGoftMode BROWNIAN $parms; }

			BROWNIAN_NO_MATSUBARA { set parms [list $NOSE_ssf::reorganization_energy($i,$j)  $NOSE_ssf::correlation_time($i,$j)];
					   NOSE_ssf::putGoftMode BROWNIAN_NO_MATSUBARA $parms; }

			BROWNIAN_LOW_TEMP_HIERARCHY { set parms [list $NOSE_ssf::reorganization_energy($i,$j)  $NOSE_ssf::correlation_time($i,$j)];
					   NOSE_ssf::putGoftMode BROWNIAN_LOW_TEMP_HIERARCHY $parms; }
			
			BROWNIAN_UNDERDAMPED { set parms [list $NOSE_ssf::reorganization_energy($i,$j)  $NOSE_ssf::brownian_omega($i,$j)];
					   NOSE_ssf::putGoftMode BROWNIAN_UNDERDAMPED $parms; }
			
			DELTA { set parms [list $NOSE_ssf::gamma($i,$j)];
					   NOSE_ssf::putGoftMode DELTA $parms; }

            GAUSSIAN { set parms [list $NOSE_ssf::reorganization_energy($i,$j)  $NOSE_ssf::correlation_time($i,$j)];
					   NOSE_ssf::putGoftMode GAUSSIAN $parms; }
		
			default {  error "Unknown CORRF TYPE"; }
		
		}
	}
	NOSE_ssf::putModeEnd
}

NOSE_ssf::putGoftEnd

# write only if the mef is present
if { $NOSE_conf::smth2source == 1 } {

	NOSE_ssf::putInitCondStart;
	NOSE_ssf::putBasisType $NOSE_ssf::basis_type(1);
	NOSE_ssf::putNrStates $NOSE_ssf::dm_nr_states(1);
	NOSE_ssf::putInitialCondition 1;
	NOSE_ssf::putInitCondEnd;

}


NOSE_ssf::putEndSections


#
# Testing module input
#

# method
set mod [NOSE_conf::getRecord module_name];
set meth [NOSE_conf::getRecord moduleMethod];

if { [string compare  $mod "QME"] == 0 } {
	
	NOSE_conf::putRecord excPosition;
	NOSE_conf::putRecord molecule;
	NOSE_conf::putRecord pathway;
	NOSE_conf::putRecord runs;
	NOSE_conf::putRecord spec_wini-dw-wst;
	NOSE_conf::putRecord localBasis;
	NOSE_conf::putRecord relaxation;
	NOSE_conf::putRecord secular;
	NOSE_conf::putRecord dephasing;
	NOSE_conf::putRecord feeding;
	NOSE_conf::putRecord draining;
	NOSE_conf::putRecord moduleMethod;
	
	if { [string compare  $meth "PT2-"] >= 0 } {
		NOSE_conf::putRecord extended_grid_extent;
	}
	
} elseif { [string compare $mod "TDPT-3"] == 0 } {
	
} elseif { [string compare $mod "MC"] == 0 } {

	NOSE_conf::putRecord moduleMethod;
	
	NOSE_conf::putRecord debug_gamma

} else {


}


NOSE_conf::putEnd


###########################################################################################
# NIS finished
###########################################################################################


#
# Start executable part of NOSE or external programs
#
if { [ string compare $dbg yes ] == 0 } {

	set result [
  	catch {
    	package require Expect;
  	} version
	];

	if { $result == 0 } {

  		spawn gdb $debug_exec;
  		interact;
  		
	} else {
		puts " "
		puts "Interactive debugging is not supported, you need to start"
		puts "debugger on the command line according to instructions below."
		puts " "
		puts "Debugging instructions: "
		puts " "
		puts "(1) Run the debugger by typing         : gdb  $debug_exec"
		puts "(2) Optionally set initial break       : break main "
		puts "(3) Run NOSE in the debugger by typing : run    "
		puts " "

	}

} else {


	if { $nis_only == yes } {
	
	
	} else {
	

		set out [open "| $loader $loader_options $driver $driver_options" r+ ];

		fileevent $out readable [list Reader $out]

		vwait done
		
	}
	
}


		
#
# END of the 'nose' script
#
