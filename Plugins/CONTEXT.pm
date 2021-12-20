=head1 LICENSE

(c) 2021-2022 Daniel Taliun

=head1 CONTACT

Daniel Taliun <daniel.taliun@mcgill.ca>
    
=cut

=head1 NAME

 CONTEXT

=head1 SYNOPSIS

 mv CONTEXT.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin CONTEXT

=head1 DESCRIPTION

 A VEP plugin that retrieves sequence context around variants.
 
=cut

package CONTEXT;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

sub feature_types {
  return ['Feature','Intergenic'];
}

sub get_header_info {
  my $self = shift;
  return {
    SEQ_5MER => '5-mer context sequence around SNV',
  }
}

sub run {
  my ($self, $tva) = @_;
  my $vf = $tva->variation_feature;

  if ($vf->{start} != $vf->{end}) {
     return {};
  }
 
  my $ref_slice = $vf->slice();
  my $seq_5mer = $ref_slice->subseq($vf->{start} - 2, $vf->{end} + 2, 1);

  return { SEQ_5MER => $seq_5mer };
}


1;
