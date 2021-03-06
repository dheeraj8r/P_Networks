$distacne_cutoff = 4.5;

%normal = (
		ALA	=>	55.7551,
		ARG	=>	93.7891,
		ASN	=>	73.4.097,
		ASP	=>	75.1507,
		CYS	=>	54.9528,
		GLN	=>	78.1301,
		GLU	=>	78.8288,
		GLY	=>	47.3129,
		HIS	=>	83.7357,
		ILE	=>	67.9452,
		LEU	=>	72.2517,
		LYS	=>	69.6096,
		MET	=>	69.2569,
		PHE	=>	93.3082,
		PRO	=>	51.331,
		SER	=>	61.3946,
		THR	=>	63.7075,
		TRP	=>	106.703,
		TYR	=>	100.719,
		VAL	=>	62.3673,
	     );

=h
Ala 55.7551
Arg 93.7891
Asn 73.4097
Asp 75.1507
Cys 54.9528
Gln 78.1301
Glu 78.8288
Gly 47.3129
His 83.7357
Ile 67.9452
Leu 72.2517
Lys 69.6096
Met 69.2569
Phe 93.3082
Pro 51.331
Ser 61.3946
Thr 63.7075
Trp 106.703
Tyr 100.719
Val 62.3673
=cut

print "Enter the location of the pdb\n";
$pdb_location = <STDIN>;
chomp $pdb_name;
print "Enter the value of I_min\n";
$cut_off = <STDIN>;
chomp $cut_off;

open(PDB,"$pdb_location") || die("Can't open the file\n");
@pdb_file = <PDB>;

$pdb_location =~ /(\w+)\.pdb/;
$pdb = $1; 
open(NODELIST,">Node_list_$pdb\_$cut_off") || die;
@node_list = ();

$i=0;
while($pdb_file[$i])
{
	if($pdb_file[$i] =~ /^ATOM\s+\d+\s+CB\s+(\w+)\s+(\w+)\s+(-?\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)/)
	{
		$residue_1 = $1;
		$chain_1 = $2;
		$res_num_1 = $3;
		$x_1 = $4;
		$y_1 = $5;
		$z_1 = $6;

		$residue_length = length $residue_1;
		if($residue_length == 4) {$residue_1 =~ s/^\w{1}//;}
	
		$pscn{"$chain_1\_$residue_1\_$res_num_1"} .= "$x_1 $y_1 $z_1\n";
		$pdb_file[$i+1] =~ /^ATOM\s+\d+\s+\w+\s+\w+\s+\w+\s+(-?\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)/;
		$res_num_2 = $1;

		if($res_num_1 eq $res_num_2)
		{
			push @node_list, "$chain_1\_$residue_1\_$res_num_1";
			print NODELIST "$chain_1\_$residue_1\_$res_num_1\n";
			$j=$i+1;
			while($pdb_file[$j] =~ /^ATOM\s+\d+\s+\w+\s+\w+\s+\w+\s+($res_num_1)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)/)
			{
				$x_2 = $2;
				$y_2 = $3;
				$z_2 = $4;
				$pscn{"$chain_1\_$residue_1\_$res_num_1"} .= "$x_2 $y_2 $z_2\n";
				$j++;
			}
			$i=$j;
		}
	}
	
	$i++;
}

=h
foreach $keyword (keys %pscn)
{
	print "$keyword: $pscn{$keyword}";
}
=cut

@adjacency = ();

$i=0;
while($node_list[$i])
{
	$node_list[$i] =~ /^\w+_(\w+)_-?\d+/;
	$residue_one = $1;
	#print $residue_one;
	@temp_one = $pscn{$node_list[$i]} =~ /(-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+)\n/g;

	$j=0;
	while($node_list[$j])
	{
		$node_list[$j] =~ /^\w+_(\w+)_-?\d+/;
		$residue_two = $1;
		#print $residue_two;
		if($i == $j)
		{
			$adjaceny[$i][$j] = 0;
		}

		if($i != $j && $j != $i+1 && $i != $j+1 && $i != $j+2 && $j != $i+2)
		{
			@temp_two = $pscn{$node_list[$j]} =~ /(-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+)\n/g;
			$count = 0;
			foreach $x (@temp_one)
			{
				$x =~ /(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)/;
				$x_one = $1;
				$y_one = $2;
				$z_one = $3;

				foreach $y (@temp_two)
				{
					$y =~ /(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)/;
					$x_two = $1;
					$y_two = $2;
					$z_two = $3;	

					$distance = sqrt( ($x_one-$x_two)*($x_one-$x_two) + ($y_one-$y_two)*($y_one-$y_two) + ($z_one-$z_two)*($z_one-$z_two) );
					if( $distance <= 4.5 )
					{
						$count++;
					}
				}
			}

			if($normal{$residue_one} != 0 && $normal{$residue_two} != 0)
			{
				$interaction = ($count*100)/(sqrt($normal{$residue_one}*$normal{$residue_two}));
			}

			if( $interaction >= $cut_off )
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























