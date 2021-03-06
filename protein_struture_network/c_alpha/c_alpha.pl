print "Enter the location of the pdb\n";
$pdb_location = <STDIN>;
chomp $pdb_name;
print "Enter C-alpha C-alpha distance cut-off\n";
$cut_off = <STDIN>;
chomp $cut_off;

open(PDB,"$pdb_location") || die("Can't open the file\n");
@pdb_file = <PDB>;
@c_alpha = ();
$i=0;
while($pdb_file[$i])
{
	if($pdb_file[$i] =~ /^ATOM\s+\d+\s+CA\s+\w+\s+\w+\s+(-?\d+)\s+/)
	{
		push @c_alpha,$pdb_file[$i];
	}

	if($pdb_file[$i] =~ /^ATOM\s+\d+\s+CA\s+\w+\s+\w+\s+(-?\d+)\s+/ && $pdb_file[$i+1] =~ /^ATOM\s+\d+\s+CA\s+\w+\s+\w+\s+(-?\d+)\s+/)
	{
		$i++;
	}
	
	$i++;
}

$pdb_location =~ /(\w+)\.pdb/;
$pdb = $1; 
open(NODELIST,">Node_list_$pdb\_$cut_off") || die;
@node_list = ();

$i=0;
@adjacency = ();

while($c_alpha[$i])
{
	$c_alpha[$i] =~ /^ATOM\s+\d+\s+CA\s+(\w+)\s+(\w+)\s+(-?\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+/;
	$res_name = $1;
	$res_chain = $2;
	$res_num = $3;
	$x_one = $4;
	$y_one = $5;
	$z_one = $6;

	push @node_list,"$res_chain $res_name $res_num";
	print NODELIST "$res_chain $res_name $res_num\n";

	$j=0;
	while($c_alpha[$j])
	{

		if($i == $j)
		{
			$adjacency[$i][$j] = 0;
		}		

		if($i != $j && $j != $i+1 && $i != $j+1)
		{
			$c_alpha[$j] =~ /^ATOM\s+\d+\s+CA\s+(\w+)\s+(\w+)\s+(-?\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+/;
			$res_name_two = $1;
			$res_chain_two = $2;
			$res_num_two = $3;
			$x_two = $4;
			$y_two = $5;
			$z_two = $6;

			$distance = sqrt( ($x_one-$x_two)*($x_one-$x_two) + ($y_one-$y_two)*($y_one-$y_two) + ($z_one-$z_two)*($z_one-$z_two) );
			if($distance <= $cut_off)
			{
				$adjacency[$i][$j] = 1;
			}

			else
			{
				$adjacency[$i][$j] = 0;
			}
		}

		else
		{
			$adjacency[$i][$j] = 0;
		}

		$j++;
	}

	$i++;
}

close NODELIST;

open(ADJACENCY,">$pdb.$cut_off.ajdm");
foreach my $x (@adjacency) 
{
	$row_string = "";
 	foreach my $y (@$x)
	{
		$row_string .= "$y "; 
	}

	print ADJACENCY "$row_string\n";
}
close ADJACENCY;

#####################################hubs###################################
@hubs = &hub(4);

sub hub
{
	my ($min_degree) = @_;
	my ($i) = 0;
	my ($j) = 0;
	my (@hubs) = ();
	my (@row_elements) = ();
	my $count;
	foreach my $row (@adjacency)
	{
		my $string = "";
		$string .= "$node_list[$i]\t";

		$j = 0;
		$count = 0;
		foreach my $row_element (@$row)
		{
			if($row_element == 1)
			{
				$string .= "$node_list[$j]\t";
			}
			$j++;
		}

		$string =~ s/\t$//;
		@row_elements = split(/\t/,$string);

		if( $min_degree <= $#row_elements )
		{
			$string =~ s/^\w+\s+\w+\s+-?\d+\t//;
			$string1 = "$#row_elements $string\n";
			push @hubs,$string1;
		}

		$i++;
	}

	open(HUBS,">$pdb.$cut_off.hubs") || die;
	print HUBS @hubs;
	close HUBS;
	return(@hubs);

}

#############################Cliques and Communities##########################
&clique;

sub clique
{
	my ($i) = 0;
	my ($j) = 0;
	my (@cfinder_input) = ();
	my $count;
	my $string1;
	foreach my $row (@adjacency)
	{
		my $string = "";
		$node_list[$i] =~ /^(\w+)\s+(\w+)\s+(-?\d+)/;
		$string = "$1_$2_$3";

		$j = 0;
		foreach my $row_element (@$row)
		{
			
			if($row_element == 1 && $j>$i)
			{
				$node_list[$j] =~ /^(\w+)\s+(\w+)\s+(-?\d+)/;
				$string1 = "$string $1_$2_$3\n";
				push @cfinder_input,$string1;
			}
			$j++;
		}
		$i++;

	}

	open(CFINDER,">$pdb.$cut_off.cfi") || die;
	print CFINDER @cfinder_input;
	close CFINDER;
	my $command_line = "./CFinder -i $pdb.$cut_off.cfi";
	`$command_line1`;
	`$command_line`;
}

##############################Cluster-Depth First Search Algorithm##############################
@marks = ();
$i=0;
while($node_list[$i])
{
	$marks[$i] = 0;
	$i++;
}

$components = 0;
for($i=0;$i<=$#node_list;$i++)
{
	if($marks[$i] == 0)
	{
		$components++;
		&dfs($i,$components);
	}

}

for($k=0;$k<=$#marks;$k++)
{
	$hash_marks{$marks[$k]} .= "$node_list[$k] ";
}

open(CLUSTER,">$pdb.$cut_off.cluster");
foreach $capword  (sort{$a <=> $b} keys(%hash_marks))
{
	print CLUSTER "$capword: $hash_marks{$capword}\n";
}
close CLUSTER;


sub dfs
{
	my ($vertex,$components) = @_;
	$marks[$vertex] = $components;

	for(my $z=0;$z<=$#node_list;$z++)
	{
		if($adjacency[$vertex][$z] == 1)
		{
			if($marks[$z] == 0)
			{
				&dfs($z,$components);
			}
		}
	}
	return;
}

