#!/usr/bin/perl -w
use 5.026;
use strict;

# invert matrix
# An Efficient and Simple Algorithm for Matrix Inversion
# Ahmad Farooq, King Khalid University, Saudi Arabia
# Khan Hamid, National University of Computer and Emerging Sciences (NUCES),
# Pakistan

# example
my $m=[[2,5],[1,3]];                 # matrix to invert
print join(', ', invert($m)), "\n";  # invert $m, printing result

sub invert($)
{
    my $m = shift;          # matrix is an array of rows
    my ($pp, $det);
    my ($rp, $pe);
    my $n = scalar(@$m);

    for ($pp = 0, $det = 1.0; $pp < $n; ++$pp) {
        $rp = $m->[$pp];        # pivot row
        $pe = $rp->[$pp];       # pivot element
        print "pe[$pp] == $pe\n";
        last if ($pe == 0);      # Epsilon test?

        $det *= $pe;
        # calculate pivot row
        for (my $j = 0; $j < $n; ++$j) {
            next if ($j == $pp);

            $rp->[$j] /= $pe;
        }

        pm($m, "pivot row $pp");
        # calculate pivot column
        for (my $i = 0; $i < $n; ++$i) {
            next if ($i == $pp);

            $m->[$i]->[$pp] /= -$pe;
        }

        pm($m, "pivot column $pp");
        for (my $j = 0; $j < $n; ++$j) {
            next if ($j == $pp);

            for (my ($i, $rj) = (0, $m->[$j]); $i < $n; ++$i) {
                next if ($i == $pp);

                $rj->[$i] += $rp->[$j] * $m->[$i]->[$pp];
            }
        }

        pm($m, "rest $pp");
        $rp->[$pp] = 1.0 / $pe;
        pm($m, "pivot $pp");
    }

    return ($pe != 0.0, $det);
}

# print matrix
sub pm($;$)
{
    my ($m, $label) = @_;
    my $n = scalar(@$m);
    print "($label) " if ($label);
    print "${n}x${n}:\n";
    for (my $i = 0; $i < $n; ++$i) {
        for (my $j = 0; $j < $n; ++$j) {
            if (defined(my $v = $m->[$i]->[$j])) {
                printf('%8.3f', $v);
            } else {
                print ' ???????';
            }
        }

        print "\n";
    }
}

