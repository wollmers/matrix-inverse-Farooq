#!perl
#use 5.026;
use strict;
use warnings;

use Math::Matrix;

use Data::Dumper;

use Test::More;



my $VERBOSE = 1;

# invert matrix
# An Efficient and Simple Algorithm for Matrix Inversion
# Ahmad Farooq, King Khalid University, Saudi Arabia
# Khan Hamid, National University of Computer and Emerging Sciences (NUCES),
# Pakistan


# A ... square matrix              $m
# a ... element                    $m->[]
# n ... dimension                  $n
# d ... determinant of the matrix  $det
sub invert($)
{
    my $A = shift;          # matrix is an array of rows
    return(undef,undef,undef)
      if ( scalar(@$A) != scalar(@{$A->[0]}) );

    my $n = scalar(@$A);

    pm($A, 'Step 0: start') if $VERBOSE;
    # Step 1: Let p = 0, d = 1;
    my $p   = 0;
    my $det = 1.0;
    # Step 2: p <= p +1
    for (my $pi = 0,$det = 1.0; $pi < $n; ++$pi) {
        $p = $pi;
        print 'Step 2: ','$p: ',$p,' $det: ',$det, ' $A->[$p]->[$p] == ',$A->[$p]->[$p],"\n" if $VERBOSE > 1;

        # Step 3: If a[p,p] == 0 then cannot calculate inverse, go to step 10.
        if ($A->[$p]->[$p] == 0) {
            print 'Step 3 ($A->[$p]->[$p] == 0): ','$A->[$p]->[$p]: ',$A->[$p]->[$p],
              ' $p: ',$p,"\n" if $VERBOSE > 1;
            last;
        }

        # Step 4: d <= d x a[p, p]
        $det = $det * $A->[$p]->[$p];
        print 'Step 4: ','$det: ',$det,"\n" if $VERBOSE > 1;

        # calculate pivot row
        # Step 5: Calculate the new elements of the pivot row by:
        #   a_new[p,j] <= a[p,j] / a[p,p] where  j = 1 .. n, j != p
        STEP5: for (my $j = 0; $j < $n; ++$j) {
            print 'Step 5: ','$j: ',$j,' $p: ',$p,"\n" if $VERBOSE > 1;
            # next if ($j == $p);
            if ($j == $p) {
              print 'Step 5 ($j == $p): ','$j: ',$j,' $p: ',$p,"\n" if $VERBOSE > 1;
              next STEP5;
            }

            $A->[$p]->[$j] = $A->[$p]->[$j] / $A->[$p]->[$p];
        }
        pm($A, 'Step 5: pivot row $p: '.$p) if $VERBOSE > 1;

        # calculate pivot column
        # Step 6: Calculate the new elements of the pivot column by:
        #   a_new[i,p] <= -(a[i,p] / a[p,p]) where  i = 1 .. n, i != p
        STEP6: for (my $i = 0; $i < $n; ++$i) {
            print 'Step 6: ','$i: ',$i,' $p: ',$p,"\n" if $VERBOSE > 1;

            if ($i == $p) {
              print 'Step 6 ($i == $p): ','$i: ',$i,' $p: ',$p,"\n" if $VERBOSE > 1;
              next STEP6;
            }

            $A->[$i]->[$p] = -($A->[$i]->[$p] / $A->[$p]->[$p]);

        }
        pm($A, 'Step 6: pivot column $p: '.$p) if $VERBOSE > 1;

        # Step 7: Calculate the rest of the new elements by:
        #   a_new[i,j] <= a[i,j] + a[p,j] x a_new[i,p]
        #     where i = 1 .. n, j = 1 .. n, & i,j != p

        OUTER7: for (my $i = 0; $i < $n; ++$i) {
            if ($i == $p) {
                print 'Step 7 ($i == $p): ','$i: ',$i,' $p: ',$p,"\n" if $VERBOSE > 1;
                next OUTER7;
            }
            INNER7: for (my $j = 0; $j < $n; ++$j) {
                if ($j == $p) {
                    print 'Step 7 ($j == $p): ','$j: ',$j,' $p: ',$p,"\n" if $VERBOSE > 1;
                    next INNER7;
                }
                # Note that in step 7 of the following algorithm a'[i, p]
                # on the LHS means that the latest value of the pivot row
                # is to be used in the calculations.
                $A->[$i]->[$j] = $A->[$i]->[$j] + $A->[$p]->[$j] * $A->[$i]->[$p];
            }
        }
        pm($A, 'Step 7: rest $p: '.$p) if $VERBOSE > 1;

        # Step 8: Calculate the new value of the current pivot location:
        #   a_new[p,p] <= 1 / a_new[p,p]
        $A->[$p]->[$p] = 1.0 / $A->[$p]->[$p];
        pm($A, 'Step 8: pivot $p: '.$p) if $VERBOSE > 1;

        # Step 9: If p < n go to step 2 (n the dimension of the matrix A).
    }

    # Step 10: Stop. If inverse exists, A contains the inverse and d is the determinant.
    pm($A, 'Step 10: stop') if $VERBOSE > 1;
    print 'Step 10: ',' $p: ',$p,' $det: ',$det,"\n" if $VERBOSE > 1;

    if ($A->[$p]->[$p] != 0.0) {
        #return ($A->[$p]->[$p] != 0.0, $det);
        return ($A->[$p]->[$p] != 0.0, $det, $A);
    }
    return ($A->[$p]->[$p] != 0.0);
}


sub invert_corr_debugging($)
{
    my $A = shift;          # matrix is an array of rows
    return(undef,undef,undef)
      if ( scalar(@$A) != scalar(@{$A->[0]}) );

    my $n = scalar(@$A);

    pm($A, 'Step 0: start') if $VERBOSE > 1;
    # Step 1: Let p = 0, d = 1;
    my $p   = 0;
    my $det = 1.0;
    # Step 2: p <= p +1
    for (my $pi = 0,$det = 1.0; $pi < $n; ++$pi) {
        $p = $pi;
        print 'Step 2: ','$p: ',$p,' $det: ',$det, ' $A->[$p]->[$p] == ',$A->[$p]->[$p],"\n" if $VERBOSE > 1;

        # Step 3: If a[p,p] == 0 then cannot calculate inverse, go to step 10.
        if ($A->[$p]->[$p] == 0) {
            print 'Step 3 ($A->[$p]->[$p] == 0): ','$A->[$p]->[$p]: ',$A->[$p]->[$p],
              ' $p: ',$p,"\n" if $VERBOSE > 1;
            last;
        }

        # Step 4: d <= d x a[p, p]
        $det = $det * $A->[$p]->[$p];
        print 'Step 4: ','$det: ',$det,"\n" if $VERBOSE > 1;



        # calculate pivot column
        # Step 6: Calculate the new elements of the pivot column by:
        #   a_new[i,p] <= -(a[i,p] / a[p,p]) where  i = 1 .. n, i != p
        STEP6: for (my $i = 0; $i < $n; ++$i) {
            print 'Step 6: ','$i: ',$i,' $p: ',$p,"\n" if $VERBOSE > 1;

            if ($i == $p) {
              print 'Step 6 ($i == $p): ','$i: ',$i,' $p: ',$p,"\n" if $VERBOSE > 1;
              next STEP6;
            }

            $A->[$i]->[$p] = -($A->[$i]->[$p] / $A->[$p]->[$p]);

        }
        pm($A, 'Step 6: pivot column $p: '.$p) if $VERBOSE > 1;

        # Step 7: Calculate the rest of the new elements by:
        #   a_new[i,j] <= a[i,j] + a[p,j] x a_new[i,p]
        #     where i = 1 .. n, j = 1 .. n, & i,j != p

        OUTER7: for (my $i = 0; $i < $n; ++$i) {
            if ($i == $p) {
                print 'Step 7 ($i == $p): ','$i: ',$i,' $p: ',$p,"\n" if $VERBOSE > 1;
                next OUTER7;
            }
            INNER7: for (my $j = 0; $j < $n; ++$j) {
                if ($j == $p) {
                    print 'Step 7 ($j == $p): ','$j: ',$j,' $p: ',$p,"\n" if $VERBOSE > 1;
                    next INNER7;
                }
                # Note that in step 7 of the following algorithm a'[i, p]
                # on the LHS means that the latest value of the pivot row
                # is to be used in the calculations.
                $A->[$i]->[$j] = $A->[$i]->[$j] + $A->[$p]->[$j] * $A->[$i]->[$p];
            }
        }
        pm($A, 'Step 7: rest $p: '.$p) if $VERBOSE > 1;

        # calculate pivot row
        # Step 5: Calculate the new elements of the pivot row by:
        #   a_new[p,j] <= a[p,j] / a[p,p] where  j = 1 .. n, j != p
        STEP5: for (my $j = 0; $j < $n; ++$j) {
            print 'Step 5: ','$j: ',$j,' $p: ',$p,"\n" if $VERBOSE > 1;
            # next if ($j == $p);
            if ($j == $p) {
              print 'Step 5 ($j == $p): ','$j: ',$j,' $p: ',$p,"\n" if $VERBOSE > 1;
              next STEP5;
            }

            $A->[$p]->[$j] = $A->[$p]->[$j] / $A->[$p]->[$p];
        }
        pm($A, 'Step 5: pivot row $p: '.$p) if $VERBOSE > 1;

        # Step 8: Calculate the new value of the current pivot location:
        #   a_new[p,p] <= 1 / a_new[p,p]
        $A->[$p]->[$p] = 1.0 / $A->[$p]->[$p];
        pm($A, 'Step 8: pivot $p: '.$p) if $VERBOSE > 1;

        # Step 9: If p < n go to step 2 (n the dimension of the matrix A).
    }

    # Step 10: Stop. If inverse exists, A contains the inverse and d is the determinant.
    pm($A, 'Step 10: stop') if $VERBOSE > 1;
    print 'Step 10: ',' $p: ',$p,' $det: ',$det,"\n" if $VERBOSE > 1;

    if ($A->[$p]->[$p] != 0.0) {
        #return ($A->[$p]->[$p] != 0.0, $det);
        return ($A->[$p]->[$p] != 0.0, $det, $A);
    }
    return ($A->[$p]->[$p] != 0.0);
}

sub invert_corr($) {
    my $A = shift;          # matrix is an array of rows
    my $n = scalar(@$A);

    # Step 1: Let p = 0, d = 1;
    my $p   = 0;
    my $det;
    # Step 2: p <= p +1
    for (my $pi = 0,$det = 1.0; $pi < $n; ++$pi) {
        $p = $pi;

        # Step 3: If a[p,p] == 0 then cannot calculate inverse, go to step 10.
        if ($A->[$p]->[$p] == 0) { last; }

        # Step 4: d <= d x a[p, p]
        $det = $det * $A->[$p]->[$p];

        # Step 6: Calculate the new elements of the pivot column by:
        #   a_new[i,p] <= -(a[i,p] / a[p,p]) where  i = 1 .. n, i != p
        STEP6: for (my $i = 0; $i < $n; ++$i) {
            if ($i == $p) { next STEP6; }
            $A->[$i]->[$p] = -($A->[$i]->[$p] / $A->[$p]->[$p]);
        }

        # Step 7: Calculate the rest of the new elements by:
        #   a_new[i,j] <= a[i,j] + a[p,j] x a_new[i,p]
        #     where i = 1 .. n, j = 1 .. n, & i,j != p
        OUTER7: for (my $i = 0; $i < $n; ++$i) {
            if ($i == $p) { next OUTER7; }
            INNER7: for (my $j = 0; $j < $n; ++$j) {
                if ($j == $p) { next INNER7; }
                # Note that in step 7 of the following algorithm a'[i, p]
                # on the LHS means that the latest value of the pivot row
                # is to be used in the calculations.
                $A->[$i]->[$j] = $A->[$i]->[$j] + $A->[$p]->[$j] * $A->[$i]->[$p];
            }
        }

        # Step 5: Calculate the new elements of the pivot row by:
        #   a_new[p,j] <= a[p,j] / a[p,p] where  j = 1 .. n, j != p
        STEP5: for (my $j = 0; $j < $n; ++$j) {
            # next if ($j == $p);
            if ($j == $p) { next STEP5; }
            $A->[$p]->[$j] = $A->[$p]->[$j] / $A->[$p]->[$p];
        }

        # Step 8: Calculate the new value of the current pivot location:
        #   a_new[p,p] <= 1 / a_new[p,p]
        $A->[$p]->[$p] = 1.0 / $A->[$p]->[$p];

        # Step 9: If p < n go to step 2 (n the dimension of the matrix A).
    }

    # Step 10: Stop. If inverse exists, A contains the inverse and d is the determinant.
    if ($A->[$p]->[$p] != 0.0) {
        return ($A->[$p]->[$p] != 0.0, $det, $A);
    }
    return ($A->[$p]->[$p] != 0.0);
}


# by Håkon Hægland
# https://stackoverflow.com/questions/57666611/bad-matrix-inversion-algorithm-or-implemented-incorrectly
sub invert_hakon($)
{
    my $m = shift;          # matrix is an array of rows
    my ($pp, $det);
    my ($rp, $pe);
    my $n = scalar(@$m);

    for ($pp = 0, $det = 1.0; $pp < $n; ++$pp) {
        $rp = $m->[$pp];        # pivot row
        $pe = $rp->[$pp];       # pivot element
        last if ($pe == 0);      # Epsilon test?

        $det *= $pe;
        # calculate pivot column
        for (my $i = 0; $i < $n; ++$i) {
            next if ($i == $pp);
            $m->[$i][$pp] /= -$pe;
        }
        for (my $j = 0; $j < $n; ++$j) { # row index
            next if ($j == $pp);
            for (my ($i, $rj) = (0, $m->[$j]); $i < $n; ++$i) {
                next if ($i == $pp);
                $rj->[$i] += $rp->[$i] * $m->[$j]->[$pp];
            }
        }
        # calculate pivot row
        for (my $j = 0; $j < $n; ++$j) {
            next if ($j == $pp);
            $rp->[$j] /= $pe;
        }
        $rp->[$pp] = 1.0 / $pe;
    }

    return ($pe != 0.0, $det, $m);
}

# print matrix
sub pm($;$)
{
    my ($A, $label) = @_;
    my $n = scalar(@$A);
    print "($label) " if ($label);
    print "${n}x${n}:\n";
    for (my $i = 0; $i < $n; ++$i) {
        for (my $j = 0; $j < $n; ++$j) {
            if (defined(my $v = $A->[$i]->[$j])) {
                printf('%8.3f', $v);
            } else {
                print ' undef';
            }
        }

        print "\n";
    }
}

sub matmult {
  my ($A,$B) = @_;                        # array of rows
  my $l = scalar( @$A );                     # l ... rows of C, i, rows A
  my $n = scalar( @{$B->[0]} );              # n ... cols of C, k, cols B
  my $m = scalar( @$B );                     # m ... j, cols A, rows B
  my $C = [];              # zero C
  for (my $i = 0;$i < $l;$i++) {          # rows of C
    for (my $k = 0;$k < $n;$k++) {        # cols of  C
      for (my $j = 0;$j < $m;$j++) {      # cols A / rows B
        $C->[$i]->[$k] += $A->[$i]->[$j] * $B->[$j]->[$k] # productsum
      }
    }
  }
  return $C;
}

sub matclone {
  my ($A) = @_;

  my $C = [];
  my $l = scalar( @$A );

  for (my $i = 0;$i < $l;$i++) {
      for (my $j = 0;$j < $l;$j++) {
        $C->[$i]->[$j] = $A->[$i]->[$j];
      }
  }
  return $C;
}



my $examples = {
  '01_wiki' => {
    'A' => [
	  [ 2, 5],
	  [ 1, 3]
    ],
    'Ainv' => [
	  [ 3,-5],
	  [-1, 2]
    ],
    'I' => [
	  [ 1, 0],
	  [ 0, 1]
    ],
    'det' => 1,
  },
  '02_wiki' => {
    'A' => [
	  [ 2, 3],
	  [ 1, 2]
    ],
    'Ainv' => [
	  [ 2,-3],
	  [-1, 2]
    ],
    'I' => [
	  [ 1, 0],
	  [ 0, 1]
    ],
    'det' => 1,
  },
  '03_author_1' => {
    'A' => [
	  [ 1, 1, 3],
	  [ 1, 3,-3],
	  [-2,-4,-4],
    ],
    'Ainv' => [
	  [    3,    1,  3/2],
	  [ -5/4, -1/4, -3/4],
	  [ -1/4, -1/4, -1/4],
    ],
    'I' => [
	  [ 1, 0, 0],
	  [ 0, 1, 0],
	  [ 0, 0, 1],
    ],
    'det' => -8,
  },

};

if (1) {
# sanity of test examples
for my $example (sort keys %$examples) {
   my $A    = $examples->{$example}->{'A'};
   my $Ainv = $examples->{$example}->{'Ainv'};
   my $I    = $examples->{$example}->{'I'};
   my $C = matmult($A, $Ainv);

   is_deeply $C, $I, 'matmult '.$example
         or diag pm($C,"got: "), pm($I,"expected: ");
}
}

if (0) {
my $max_tests   = 1;
my $count_tests = 0;

for my $example (sort keys %$examples) {
   $count_tests++;
   last if ($count_tests > $max_tests);

   my $A    = $examples->{$example}->{'A'};
   my $Ainv = $examples->{$example}->{'Ainv'};
   my $I    = $examples->{$example}->{'I'};

   print '*** invert '.$example,"\n" if $VERBOSE;

   my ($p_ele_ok, $det, $C) = invert($A);

   print '$C: ',pm($C,'invert '.$example) if $VERBOSE;

   is_deeply $C, $Ainv, 'invert '.$example
        or diag pm($C,"got: "), pm($Ainv,"expected: ");

   is $det, $examples->{$example}->{'det'}, 'invert det: '.$example.'$det';
}
}

if (1) {
my $max_tests   = 99;
my $count_tests = 0;

for my $example (sort keys %$examples) {
   $count_tests++;
   last if ($count_tests > $max_tests);

   my $A    = matclone($examples->{$example}->{'A'});
   my $Ainv = matclone($examples->{$example}->{'Ainv'});
   my $I    = matclone($examples->{$example}->{'I'});

   print '*** invert_corr '.$example,"\n" if $VERBOSE > 1;

   pm($A,'invert_corr '.$example.' input $A') if $VERBOSE;

   my ($p_ele_ok, $det, $C) = invert_corr($A);

   pm($C,'invert_corr '.$example.' result $C') if $VERBOSE;

   is_deeply $C, $Ainv, $example.' invert_corr Ainv '
        or diag pm($C,"got: "), pm($Ainv,"expected: ");

   is $det, $examples->{$example}->{'det'}, $example.' invert_corr'.' det: '.$det;
}
}

if (1) {
my $max_tests   = 99;
my $count_tests = 0;

for my $example (sort keys %$examples) {
   $count_tests++;
   last if ($count_tests > $max_tests);

   my $A    = matclone($examples->{$example}->{'A'});
   my $Ainv = matclone($examples->{$example}->{'Ainv'});
   my $I    = matclone($examples->{$example}->{'I'});

   print '*** invert_hakon '.$example,"\n" if $VERBOSE > 1;

   pm($A,'invert_hakon '.$example.' input $A') if $VERBOSE;

   my ($p_ele_ok, $det, $C) = invert_hakon($A);

   pm($C,'invert_hakon '.$example.' result $C') if $VERBOSE;

   is_deeply $C, $Ainv, $example.' invert_hakon Ainv'
        or diag pm($C,"got: "), pm($Ainv,"expected: ");

   is $det, $examples->{$example}->{'det'}, $example.' invert_hakon'.' det: '.$det;
}
}

if (1) {
# test examples with Math::Matrix as a reference
for my $example (sort keys %$examples) {
   my $A    = matclone($examples->{$example}->{'A'});
   my $Ainv = matclone($examples->{$example}->{'Ainv'});
   my $I    = matclone($examples->{$example}->{'I'});

   print '*** Math::Matrix->invert '.$example,"\n" if $VERBOSE > 1;

   pm($A,'Math::Matrix->invert '.$example.' input $A') if $VERBOSE;

   my $A_obj = Math::Matrix->new(@$A);

   my $C    = [@{$A_obj->invert}]; # $A_obj->invert returns a blessed object

   pm($C,'Math::Matrix->invert '.$example.' result $C') if $VERBOSE;

   is_deeply $C, $Ainv, $example.' Math::Matrix->invert Ainv'
        or diag pm($C,"got: "), pm($Ainv,"expected: ");

   is $A_obj->determinant, $examples->{$example}->{'det'},
       $example.' Math::Matrix->det: '.$examples->{$example}->{'det'};
}
}

done_testing;
