#!/usr/local/bin/perl
package Statistics::TTest;
use strict;
use Carp;
use POSIX;
use vars qw($VERSION $AUTOLOAD);
use Statistics::Distributions qw (tdistr fdistr tprob fprob);
use Statistics::PointEstimation;
$VERSION='1.1';
my %fields=
(
	's1'          =>undef,     #sample 1 
	's2'          =>undef,     #sample 2 
	'significance' =>undef,
	'alpha' =>undef,
	'mean_difference' =>undef,
	'variance' =>undef,
	'standard_deviation' =>undef,
	'standard_error' =>undef,
	'standard_error_equal' =>undef,
	'standard_error_unequal' =>undef,
	'f_cutoff' =>undef,
	'f_statistic' =>undef,
	'df1' =>undef,
	'df2' =>undef,
	'df' => undef,
	'df_equal' =>undef,
	'df_unequal' =>undef,
	'equal_variance'=>undef,
	't_value' =>undef,
	't_statistic' =>undef,
	't_statistic_equal' =>undef,
	't_statistic_unequal' =>undef,
	't_prob' =>undef,
	'null_hypothesis' =>undef,
	'delta' =>undef,
	'lower_clm' =>undef,
	'upper_clm' =>undef,
	'valid' =>undef
);


sub new
{
	my $proto = shift;
        my $class = ref($proto) || $proto;
	my $self= {%fields};
	bless($self,$class);
	return $self;
}

sub load_data
{
	my $self=shift;
	undef $self->{valid};
	my ($sample1,$sample2)=@_;
	if((ref $sample1 ne 'ARRAY')||(ref $sample2 ne 'ARRAY'))
	{
		croak "Invalid input type for load_data() function.\n the 2 inputs must be array references.";
	}

	my $s1= new Statistics::PointEstimation;
	my $s2= new Statistics::PointEstimation;

	if($self->significance())
	{
		$s1->set_significance($self->significance());
		$s2->set_significance($self->significance());
	}

	$s1->add_data($sample1);
	$s2->add_data($sample2);

	croak "Invalid sample size.  sample size for the 2 samples are",$s1->count," and ",$s2->count() ,
	      ". the sample size must be greater than 1" if(($s1->count() <=1)||($s2->count() <=1 ));
	$self->{'s1'}=$s1;
	$self->{'s2'}=$s2;

	return $self->perform_t_test();
	


}

sub perform_t_test
{
	my $self=shift;
	my $s1=$self->s1();
	my $s2=$self->s2();
	$self->{valid}=0;
	$self->{significance}=95 if (!defined($self->{significance}));
	$self->{alpha}=(100-$self->{significance})/2;
        $self->{alpha}/=100;
	$self->{mean_difference}=$s1->{mean}-$s2->{mean};
	$self->{variance}=($s1->{df}*$s1->{variance}+$s2->{df}*$s2->{variance})/($s1->{df}+$s2->{df});
	$self->{standard_deviation}=sqrt($self->{variance});
	$self->{standard_error_equal}=sqrt(1/$s1->{count}+1/$s2->{count})*$self->{standard_deviation};
	$self->{standard_error_unequal}=sqrt(($s1->{standard_error})**2+($s2->{standard_error})**2);

	$self->{df_equal}=$s1->{df}+$s2->{df};
	$self->{df_unequal}= (($s1->{standard_error})**4/$s1->{df}+($s2->{standard_error})**4/$s2->{df})?
				 ((($s1->{standard_error})**2+($s2->{standard_error})**2)**2
				   /
			           (($s1->{standard_error})**4/$s1->{df}+($s2->{standard_error})**4/$s2->{df})) :
				 ($self->{df_equal});

	$self->{t_statistic_equal}=($self->{standard_error_equal})?
				   (abs ($self->{mean_difference}/$self->{standard_error_equal})):99;
	$self->{t_statistic_unequal}=($self->{standard_error_unequal})?
				   (abs ($self->{mean_difference}/$self->{standard_error_unequal})):99;
	my ($df1,$df2);
	if($s1->{variance}>=$s2->{variance})
	{
		$df1=$s1->{df};
		$df2=$s2->{df};
		$self->{f_statistic}=($s2->{variance})?($s1->{variance}/$s2->{variance}):99;

	}
	else
	{
	
		$df1=$s2->{df};
		$df2=$s1->{df};
		$self->{f_statistic}=($s1->{variance})?($s2->{variance}/$s1->{variance}):99;
	}
	($self->{df1},$self->{df2})=($df1,$df2);
	$self->{f_cutoff}=fdistr($df1,$df2,$self->{alpha});

	if($self->{f_statistic}<=$self->{f_cutoff})
	{ $self->{equal_variance}=1; }
	else
	{ $self->{equal_variance}=0; }


	if($self->{equal_variance})
	{
		$self->{standard_error}=$self->{standard_error_equal};
		$self->{t_statistic}=$self->{t_statistic_equal};
		$self->{df}=$self->{df_equal};
                                                                                                                                                                                                                                                                                                                                      
	}
	else
	{
		$self->{standard_error}=$self->{standard_error_unequal};
		$self->{t_statistic}=$self->{t_statistic_unequal};
		$self->{df}=$self->{df_unequal};

	}
	$self->{t_prob}=1- abs (tprob(floor($self->{df}),-1*$self->{t_statistic})-tprob(floor($self->{df}),$self->{t_statistic}));
	if($self->{t_prob}<$self->{alpha}*2)
	{
		$self->{null_hypothesis}='rejected';
	}	
	else
	{
		$self->{null_hypothesis}='not rejected';
	}

	$self->{t_value}=tdistr(floor($self->df()),$self->alpha());
	$self->{delta}=$self->t_value()*$self->standard_error();
	$self->{lower_clm}=$self->{mean_difference}-$self->{delta};
	$self->{upper_clm}=$self->{mean_difference}+$self->{delta};
	$self->{valid}=1 if (($s1->{variance})&&($s2->{variance}));
	return !($self->{valid});

}

sub print_t_test
{
	my $self=shift;
	foreach my $k ( keys %$self)
        {
		next if ($k eq 's1') || ($k eq 's2'); 
                print "$k: $self->{$k} \n";
        }
        return 1;

}
sub output_t_test
{
	my $self=shift;
	my $s1=$self->{s1};
	my $s2=$self->{s2};
	croak "no data. s1 or s2 is empty or the variance =0.\n" if ((!defined($s1))||(!defined($s2))||($self->valid!=1));
	print "*****************************************************\n\n";
	$s1->output_confidence_interval('1');
	print "*****************************************************\n\n";
	$s2->output_confidence_interval('2');
	print "*****************************************************\n\n";

        print "Comparison of these 2 independent samples.\n";
        print "\t F-statistic=",$self->f_statistic()," , cutoff F-statistic=",$self->f_cutoff(),
		" with alpha level=",$self->alpha*2," and  df =(",$self->df1,",",$self->df2,")\n"; 
	if($self->{equal_variance})
	{ print "\tequal variance assumption is accepted(not rejected) since F-statistic < cutoff F-statistic\n";}
	else
	{ print "\tequal variance assumption is  rejected since F-statistic > cutoff F-statistic\n";}

	print "\tdegree of freedom=",$self->df," , t-statistic=T=",$self->t_statistic," Prob >|T|=",$self->{t_prob},"\n";
	print "\tthe null hypothesis (the 2 samples have the same mean) is ",$self->null_hypothesis(),
		 " since the alpha level is ",$self->alpha()*2,"\n";
	print "\tdifference of the mean=",$self->mean_difference(),", standard error=",$self->standard_error(),"\n";
	print "\t the estimate of the difference of the mean is ", $self->mean_difference()," +/- ",$self->delta(),"\n\t",
                " or (",$self->lower_clm()," to ",$self->upper_clm," ) with ",$self->significance," % of confidence\n"; 

	return 1;
}
sub set_significance
{
	my $self=shift;
        my $significance=shift;
        croak "Invalid Significance. $significance should be 0-100 usually 90,95,99\n"
		 unless (($significance>0)&&($significance<100));
	$self->{significance}=$significance;
        if($self->{s1}&&$self->{s2})
	{
		$self->{s1}->set_significance($significance);
		$self->{s2}->set_significance($significance);
		$self->perform_t_test();

	}
        return 1;
	
}

sub AUTOLOAD
{
	my $self = shift;
  	my $type = ref($self)
    	or croak "$self is not an object";
  	my $name = $AUTOLOAD;
  	$name =~ s/.*://;     ##Strip fully qualified-package portion
  	return if $name eq "DESTROY";
  	unless (exists $self->{$name} )
	{
    		croak "Can't access `$name' field in class $type";
  	}
	 ##Read only method 
	 return $self->{$name};
}

1;

#perform t-test using sufficient statistics
package Statistics::TTest::Sufficient;
use strict;
use Carp;
use vars qw($VERSION @ISA $AUTOLOAD);
use Statistics::PointEstimation;
use Statistics::TTest;
use POSIX;
@ISA= qw (Statistics::TTest);
$VERSION = '1.1';


sub load_data{
	my $self=shift;
        undef $self->{valid};
        my ($sample1,$sample2)=@_;

        
        my $s1= new Statistics::PointEstimation::Sufficient;
        my $s2= new Statistics::PointEstimation::Sufficient;

        if($self->significance())
        {
                $s1->set_significance($self->significance());
                $s2->set_significance($self->significance());
        }
        
        $s1->load_data($sample1->{count},$sample1->{mean},$sample1->{variance});
        $s2->load_data($sample2->{count},$sample2->{mean},$sample2->{variance});
        
        croak "Invalid sample size.  sample size for the 2 samples are",$s1->count," and ",$s2->count() ,
              ". the sample size must be greater than 1" if(($s1->count() <=1)||($s2->count() <=1 ));
        $self->{'s1'}=$s1;
        $self->{'s2'}=$s2;
                
        return $self->perform_t_test();
        

}

1;
__END__

=head1 NAME

 Statistics::TTest - Perl module to perform T-test on 2 independent samples
 Statistics::TTest::Sufficient - Perl module to perfrom T-Test on 2 indepdent samples using sufficient statistics

=head1 SYNOPSIS
 
 #example for Statistics::TTest
 use Statistics::PointEstimation;
 use Statistics::TTest;
 my @r1=();
 my @r2=();
 my $rand;
        
  for($i=1;$i<=32;$i++) #generate a uniformly distributed sample with mean=5   
  {

          $rand=rand(10);
          push @r1,$rand;
          $rand=rand(10)-2;
          push @r2,$rand;
  }


 my $ttest = new Statistics::TTest;  
 $ttest->set_significance(90);
 $ttest->load_data(\@r1,\@r2);  
 $ttest->output_t_test();   
 $ttest->set_significance(99);
 $ttest->print_t_test();  #list out t-test related data
  
 #the following thes same as calling output_t_test() (you can check if $ttest->{valid}==1 to check if the data is valid.)
 my $s1=$ttest->{s1};  #sample 1  a Statistics::PointEstimation object
 my $s2=$ttest->{s2};  #sample 2  a Statistics::PointEstimation object
 print "*****************************************************\n\n";
 $s1->output_confidence_interval('1');
 print "*****************************************************\n\n";
 $s2->output_confidence_interval('2');
 print "*****************************************************\n\n";
 
 print "Comparison of these 2 independent samples.\n";
 print "\t F-statistic=",$ttest->f_statistic()," , cutoff F-statistic=",$ttest->f_cutoff(),
 	" with alpha level=",$ttest->alpha*2," and  df =(",$ttest->df1,",",$ttest->df2,")\n"; 
 if($ttest->{equal_variance})
 { print "\tequal variance assumption is accepted(not rejected) since F-statistic < cutoff F-statistic\n";}
 else
 { print "\tequal variance assumption is  rejected since F-statistic > cutoff F-statistic\n";}
 
 print "\tdegree of freedom=",$ttest->df," , t-statistic=T=",$ttest->t_statistic," Prob >|T|=",$ttest->{t_prob},"\n";
 print "\tthe null hypothesis (the 2 samples have the same mean) is ",$ttest->null_hypothesis(),
 	 " since the alpha level is ",$ttest->alpha()*2,"\n";
 print "\tdifference of the mean=",$ttest->mean_difference(),", standard error=",$ttest->standard_error(),"\n";
 print "\t the estimate of the difference of the mean is ", $ttest->mean_difference()," +/- ",$ttest->delta(),"\n\t",
 	" or (",$ttest->lower_clm()," to ",$ttest->upper_clm," ) with ",$ttest->significance," % of confidence\n"; 
 
 #example for Statistics::TTest::Sufficient
 use Statistics::PointEstimation;
 use Statistics::TTest;
        
 my %sample1=(
        'count' =>30,
        'mean' =>3.98,
        'variance' =>2.63
                );

 my %sample2=(
        'count'=>30,
        'mean'=>3.67,
        'variance'=>1.12
        );


 my $ttest = new Statistics::TTest::Sufficient;  
 $ttest->set_significance(90);
 $ttest->load_data(\%sample1,\%sample2);  
 $ttest->output_t_test();
 #$ttest->s1->print_confidence_interval();
 $ttest->set_significance(99);
 $ttest->output_t_test();
 #$ttest->s1->print_confidence_interval();   


=head1 DESCRIPTION

=head2 Statistics::TTest

 This is the Statistical T-Test module to compare 2 independent samples. It takes 2 array of point measures, 
 compute the confidence intervals using the PointEstimation module (which is also included in this package)
  and use the T-statistic to test the null hypothesis. If the null hypothesis is rejected, the difference 
  will be given as the lower_clm and upper_clm of the TTest object. 

=head2 Statistics::TTest::Sufficient

 This module is a subclass of Statistics::TTest. Instead of taking the real data points as the input, 
 it will compute the confidence intervals based on the sufficient statistics and the sample size inputted. 
 To use this module, you need to pass the sample size, the sample mean , and the sample variance into the load_data()
 function. The output will be exactly the same as the Statistics::TTest Module.
 


=head1 AUTHOR

Yun-Fang Juan , Yahoo! Inc.  (yunfang@yahoo-inc.com)

=head1 SEE ALSO

Statistics::Descriptive Statistics::Distributions Statistics::PointEstimation

=cut 

