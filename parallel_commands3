#!/usr/bin/perl
#v3.0

# parallel_command.pl
# Adapted from parallelRepeatMasker
#
# Example:
# $Bin/parallel_command.pl  $num_cpus  file.of.commands
# e.g. file of commands looks like: primer3_script.unformatted
no warnings 'deprecated';
select STDERR; $| = 1;  # make unbuffered - standard error
select STDOUT; $| = 1;  # make unbuffered - standard output

# sub parallel_commands (@) 
 
# Get the first command line parameters.
$num_cpus        = shift @ARGV;

print STDERR "\n";
print STDERR "num_cpus        = [$num_cpus]\n";

$command_file=shift @ARGV;
print STDERR "\n";
print STDERR "command file= [$command_file]\n";

open COMMANDS,"$command_file" or die "Can't open $command_file for reading: $!\n";

# Load up the commands to run in parallel
@commands = ();
while ($command=<COMMANDS>) {
    chomp($command);
    if( ($command =~ /^$/) || ($command =~ /^#/) ) { # ignore comments and blank lines
			       next;
			   }
  push @commands, $command;
}
close COMMANDS;

# $count= 1;
# print STDERR "\n";
# foreach $command (@commands) {
  # printf STDERR "%3d. %s\n", $count++, $command;
# }


$len = scalar @commands;
print STDERR "\n";
print STDERR "commands[0] = $commands[0]\n";
print STDERR ".\n.\n.\n";
print STDERR "commands[",$len - 1,"] = $commands[$len - 1]\n";
print STDERR "\n";



$nice_value         = 0;
$schedule_algorithm = 0;
$debug              = 0;


parallel_commands (\@commands, $num_cpus, $nice_value, $schedule_algorithm, $debug);
exit 0;



# sub parallel_commands (\@$$$$) {
sub parallel_commands () {



    # Take advantage of multi-processors if available

    # If you use this kind of processing, you should make sure jobs
    # complete by checking its output.  For example: A blast file should
    # have a "cpu time" near the bottom of the file. Another example:
    # a phd file ends with the line END_SEQUENCE.  Please be careful,
    # these detectors for completion may change over time.

    use strict;

    my $FUNCTION = "parallel_commands";

    # PARAMETERS

    # $A_unix_commands    Array of commands
    # $S_processors       Number of processors requested
    # $S_nice             Level of nice to use
    # $S_schedule_type    Scheduling algorithm to use

    my ($A_unix_commands, $S_processors, $S_nice, $S_schedule_type, $DEBUG) = @_;

    my $n_cpus;         # determine the number of processors available
    my $n_processes;    # number of processes allowed
    my %child_pid;      # store subprocess pids here
    my $finished_pid;   # finished children pids
    my @commands;       # copy of commands to be pop'd
    my $command;        # command to be forked
    my $nice_command;   # complete command with nice added in front
    my $pid;            # pid of child processes forked

    my %command_num;    # Maps child pids to 1, 2, ...
    my $command_count;  # Count of commands that have been started.
    my $date;


    # PROCESSORS PARAMETER
    # number of processor to be used

    # DETERMINE NUMBER OF PROCESSORS TO BE USED
    # schedule_type options to be added into this section
    # currently using a simple version


    $n_processes = $S_processors;


    if ($DEBUG) {
	my $node;
	chomp ($node = `hostname`);
	printf STDERR "$FUNCTION: Using %d subprocess(es) on %s\n",
	    $n_processes, $node;
    }


    # NICE PARAMETER
    # you must have privileges to use NEGATIVE nice values

    if (($S_nice < -19) || ($S_nice > 19)) {
	printf STDERR "$FUNCTION: WARNING - nice parameter out of range, using 19\n";
	$S_nice = 19;
    }


    # SCHEDULING PARAMETER
    # currently not used and ignored
    # determines scheduling algorithm used to distributing commands


    # START PROCESSES
    $command_count = 0;             # Number of commands that have started.
    %child_pid = ();                # hash of PIDs running
    @commands = @$A_unix_commands;  # make copy so we can remove items from array
    chomp @commands;                # remove line terminator if present

    while (scalar(@commands)) {

	if (scalar(keys %child_pid) < $n_processes) {

	    $command = shift @commands;
#	    $nice_command = "/bin/nice -n $S_nice $command";
#	    $nice_command =   $niceX . $command;
	    $nice_command = $command;
#            print "nice_cmd =  $nice_command\n";
	    if ($pid = fork()) {
		# I am the parent/master process
		$child_pid{$pid}++;
		$command_count++;
		$date = `date`;
#		#printf STDERR "Start command %3d %s", $command_count, $date;
		$command_num{$pid} = $command_count;
		printf "[%2d] pid %d added\n", scalar(keys %child_pid), $pid if $DEBUG;
	    }
	    else {
		# I am the child/slave process
		chomp $nice_command;
		printf "[ \$] $nice_command\n" if $DEBUG;
		exec $nice_command;
		exit 0;
	    }
	} else {
	    # we've reached the number of processes limit
	    last;
	}
    }

    # WAITING FOR PROCESSES TO COMPLETE

    while (scalar(keys %child_pid)) {

	# PROCESS COMPLETED - REMOVE FROM LIST

	$finished_pid = wait();
	delete $child_pid{$finished_pid};

        $date = `date`;
#	#printf STDERR "End   command %3d %s", $command_num{$finished_pid}, $date;
	printf "[%2d] pid %d terminated\n", scalar(keys %child_pid), $finished_pid if $DEBUG;

	# FORK MORE PROCESSES (should be same code as above)

	while (scalar(@commands)) {

	    if (scalar(keys %child_pid) < $n_processes) {

		$command = shift @commands;
#		$nice_command = "nice $command";
#	        $nice_command =   $niceX . $command;
		$nice_command = $command;

		if ($pid = fork()) {
		    # I am the parent/master process
		    $child_pid{$pid}++;

		    $command_count++;
		    $date = `date`;
	#	    #printf STDERR "Start command %3d %s", $command_count, $date;
		    $command_num{$pid} = $command_count;
		    printf "[%2d] pid %d added\n", scalar(keys %child_pid), $pid if $DEBUG;
		}
		else {
		    # I am the child/slave process
		    chomp $nice_command;
		    printf "[ \$] $nice_command\n" if $DEBUG;
		    exec $nice_command;
		    exit 0;
		}
	    } else {
		# we've reached the number of processes limit
		last;
	    }
	}
    }

    return 1;
}


##########################################################################
#  Function: n_CPUs                                                      #
##########################################################################
 
sub n_CPUs () {

    # Returns the number of CPUs found on the system.
    # If the number of processors can't be determined,
    # "1" will be returned.

    use strict;

    my ($result, $CPUs);
 
    $result = `/usr/local/bin/sysinfo 2>/dev/null`;
    ($CPUs) = ($result =~ /Number of CPUs is\s+(\d+)/);
 
    if ($CPUs =~ /\d+/) {
	return $CPUs;
    } else {
	# can't determine the number of CPUs, therefore assume 1 processor
	return 1;
    }
}

