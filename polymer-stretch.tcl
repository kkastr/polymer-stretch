set N [lindex $argv 0]
set rseed  [lindex $argv 1]
set filename [lindex $argv 2]
set vis [lindex $argv 3]


t_random seed $rseed



set boxx 400
set boxy 400
set boxz 400 

set cx [expr $boxx/8.0]
set cy [expr $boxy/8.0]
set cz [expr $boxz/8.0]
set tmax 10000
set nmax 10
set temp 1; set gamma 1.0; set gamma_equilibration 0.1



setmd box_l $boxx $boxy $boxz
setmd periodic 1 1 1
setmd time_step 0.01; setmd skin 0.4
thermostat langevin $temp $gamma


#FENE potential
set k_fene 30.0
set r_fene 1.5

#Angular Potential - Harmonic
set k_angle [expr 10.0 / $temp]
set pi 3.14159

#Shifted Lennard-Jones
set eps 1.0
set sigma 1.0
set lj_cutoff 1.12246
set lj_shift 0.25
set lj_offset 0 



inter 0 fene $k_fene $r_fene
inter 1 angle $k_angle $pi
inter 0 0 lennard-jones $eps $sigma $lj_cutoff $lj_shift $lj_offset
inter 1 0 lennard-jones $eps $sigma $lj_cutoff $lj_shift $lj_offset

set t_trans 0
set trans_flag 0
set fixed_N [expr $N/2]
set equil_time [expr 10.0 * $N]
set t_pore 1
set z_line [expr $cz - $t_pore/2]
set force [expr -5.0]
set lsens 1.5
set lcav 19.5 
set rcav 99.5
set lfilt 4.5
set rsens 1.8
set rfilt 4.5


set mirrorporez [expr $cz + 4.0 * $lcav + 2.0 * $lfilt]
set z_top [expr $cz+$mirrorporez - $t_pore/2] 




part [expr $N] pos $cx $cy $cz type 99 fix
constraint pore center [expr $cx] [expr $cy] [expr $cz] axis 0 0 1 radius 1.5 length 1 type 1
constraint pore center [expr $cx] [expr $cy] [expr 2.0 * $cz + $lcav+0.5 + $l_+0.5] axis 0 0 1 radius $rcav length $lcav type 1
constraint pore center [expr $cx] [expr $cy] [expr 2.0 * $cz + 2.0 * ($lcav + 0.5) + $lsens+0.5 + $lfilt+0.5] axis 0 0 1 radius $rfilt length $lfilt type 1
constraint pore center [expr $cx] [expr $cy] [expr 2.0 * $cz + $lcav+0.5 + $l_+0.5 + $mirrorporez] axis 0 0 1 radius $rcav length $lcav type 1
constraint pore center [expr $cx] [expr $cy] [expr 2.0 * $cz + 2.0 * ($lcav + 0.5) + $lsens+0.5 + $lfilt+0.5 + $mirrorporez] axis 0 0 1 radius $rfilt length $lfilt type 1





for { set i 0 } { $i < $N } { incr i } {
	#set x [expr $cx - $N/2 + $i]
	set x [expr $cx]
	set y [expr $cy]
	set z [expr $cz  + $i]
	part $i pos $x $y $z type 0 


	if { $i > 0 } {
		part $i bond 0 [expr $i - 1]
	}
	if { $i > 1} {
		part [expr $i - 1] bond 1 $i [expr $i - 2]
	}
}

if { $vis == 1 } {
    prepare_vmd_connection vmdout
    imd listen 100
    imd positions
}



set flag 0
set t 0	
set n 0
set fail 0
set success 0
set trans_flag 0

while {$flag == 0} {

	for { set i 0 } { $i < $N } { incr i } {
		set x [expr $cx - $N/2 + $i]
		set y [expr $cy]
		set z [expr $cz  + 100]
		part $i pos $x $y $z type 0 ext_force 0 0 0 
	}

	part $fixed_N fix

	thermostat langevin $temp $gamma_equilibration
	for {set i 0} {$i < $equil_time} {incr i} {
	    
	    integrate 100
	    imd positions
	
	}
	thermostat langevin $temp $gamma

	part $fixed_N unfix

	puts "equilibrated."



	for {set i 0} {$i < $N} {incr i} {
		part $i ext_force $force 1 1
	}

	#puts "I'm back in the first while"
	while {1} {
		set z_list {}
		set r_list {}		
		if { $n > $nmax } {
			puts "n > nmax"
			set flag 1
			break
		}

		for {set i 0} { $i < $N } {incr i} {
			set x [lindex [part $i print pos] 0]
			set y [lindex [part $i print pos] 1]
			set z [lindex [part $i print pos] 2]
			set r [expr sqrt(($x-$cx)*($x-$cx) + ($y-$cy)*($y-$cy) + ($z-$cz)*($z-$cz))]
			#puts $r
			lappend z_list $z
			lappend r_list $r
			#puts "hi, I'm getting particle positions"
		}
		set z_min [::tcl::mathfunc::min {*}$z_list]
		set z_max [::tcl::mathfunc::max {*}$z_list]

		if {$z_min < $z_line} {
			
			if {$trans_flag == 0 } {
				#puts "inside translocation if"
				set t_thread $t
				

				set rg_calc_trans [analyze rg 0 1 $N]

				set trans_flag [expr $trans_flag + 1 ]
			}
		}
		#puts $z_max
		if {$z_max < $z_line} {
			puts "zmax less than zline"
			set t_last_thread $t
			puts $t_last_thread
			set t_trans [expr $t_last_thread - $t_thread]
			set trans_time [open "data/${filename}_$N/trans_time_$N-$rseed.dat" "a"]
			puts $trans_time "$t_trans $N $t_last_thread $t_thread"
			close $trans_time
			set n [expr $n + 1.0]
			set trans_flag 0
			break
		}
		if {$t_trans != 0} {
			set t_trans 0 
		}
		integrate 100
		imd positions
		
		incr t
	}
}

