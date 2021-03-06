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
	#ATOM    139  CA  GLU A  90      -9.151  -0.537  11.459  1.00  8.51           C 
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
	

#$grep_string = "grep -P '^ATOM\\s+\\d+\\s+CA' $pdb_location";
#print $grep_string;
#@c_alpha = `$grep_string`;
#print @c_alpha;

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
	#print "$res_name $res_chain $res_num $x_one $y_one $z_one\n";

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
			#print "$res_name_two $res_chain_two $res_num_two $x_two $y_two $z_two\n";

			$distance = sqrt( ($x_one-$x_two)*($x_one-$x_two) + ($y_one-$y_two)*($y_one-$y_two) + ($z_one-$z_two)*($z_one-$z_two) );
			if($distance <= $cut_off)
			{
				$adjacency[$i][$j] = 1;
				#print "$distance $i $j $adjacency[$i][$j]\n";
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

#=h
open(ADJACENCY,">Adjacency_Matrix_$pdb\_$cut_off");
foreach my $x (@adjacency) 
{
	$row_string = "";
 	foreach my $y (@$x)
	{
		$row_string .= "$y "; 
	}

	#$row_string =~ s/ $/\n/;
	print ADJACENCY "$row_string\n";
}
close ADJACENCY;

#=h
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
		#print $string."\n";
		#print @row_elements;
		@row_elements = split(/\t/,$string);
		#pop @row_elements;
		#print "$i $#row_elements\n";
		#open(HUBS,">hubs_$pdb") || die;

		if( $min_degree <= $#row_elements )
		{
			#print HUBS "@row_elements\n";
			$string =~ s/^\w+\s+\w+\s+-?\d+\t//;
			$string1 = "$#row_elements $string\n";
			push @hubs,$string1;
		}

		$i++;
	}

	#pop @row_elements
	open(HUBS,">hubs.$pdb.$cut_off") || die;
	print HUBS @hubs;
	close HUBS;
	return(@hubs);

}


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

	open(CFINDER,">Cfinder_input.$pdb.$cut_off") || die;
	print CFINDER @cfinder_input;
	close CFINDER;
	my $command_line = "./CFinder -i Cfinder_input.$pdb.$cut_off";
	$command_line1 = "rm -f Cfinder_input.$pdb.$cut_off_files";
	`$command_line1`;
	`$command_line`;
}
#=cut

@marks = ();
$i=0;
while($node_list[$i])
{
	#print "$i $node_list[$i]\n";
	$marks[$i] = 0;
	$i++;
}

$components = 0;
for($i=0;$i<=$#node_list;$i++)
{
	#print $#node_list;
	if($marks[$i] == 0)
	{
		$components++;
		#print "$i $components\n";
		&dfs($i,$components);
	}

}

for($k=0;$k<=$#marks;$k++)
{
	$hash_marks{$marks[$k]} .= "$node_list[$k] ";
}

open(CLUSTER,">cluster.$pdb.$cut_off");
foreach $capword  (sort{$a <=> $b} keys(%hash_marks))
{
	print CLUSTER "$capword: $hash_marks{$capword}\n";
}
close CLUSTER;


sub dfs
{
	($vertex,$components) = @_;
	$marks[$vertex] = $components;
	#print "\n$vertex $components $marks[$vertex]\n";

	for($z=0;$z<=$#node_list;$z++)
	{
		if($adjacency[$vertex][$z] == 1)
		{
			#print "$vertex $z\n";
			if($marks[$z] == 0)
			{
				&dfs($z,$components);
			}
		}
	}
	return;
}

#print "@marks\n";
















