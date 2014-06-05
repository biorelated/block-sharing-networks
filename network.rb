#!/usr/bin/env ruby

require 'bio-dbla-classifier'
require 'commander/import'

program :version, '0.0.1'
program :description, 'generates a network from DBL-alpha amino acid sequence tags'
program :author, 'George Githinji email: ggithinji@kemri-wellcome.org'

command :create do |c|
  #A structure to hold the name and pspb types
  class Pspb < Struct.new(:id,:name,:pspb1,:pspb2,:pspb3,:pspb4)
    def all_pspbs
      [pspb1,pspb2,pspb3,pspb4]
    end
  end

  c.syntax = 'network create [options]'
  c.summary = 'creates a network given a list of amino acid sequences'
  c.description = 'This command generates a tab-separated list of 2 columns where each entry in each column is a node.'
  c.example 'network', 'create --infile FILENAME --blocksize INTEGER  --outfile FILENAME'
  c.option '--infile FILE', 'specify an input file(amino acid only)'
  c.option '--blocksize INTEGER', Integer, 'specify the length of the pspb block. Default is 10'

  c.action do |args, options|
    options.default :blocksize => 10
    seq_file   = options.infile
    blocksize = options.blocksize

    pspbs = Bio::FlatFile.open(seq_file).each_with_index.map do |entry,index|
      pspb1 = Bio::Sequence::AA.new(entry.seq).pspb1(0,blocksize)
      pspb2 = Bio::Sequence::AA.new(entry.seq).pspb2(0,blocksize)
      pspb3 = Bio::Sequence::AA.new(entry.seq).pspb3(0,blocksize)
      pspb4 = Bio::Sequence::AA.new(entry.seq).pspb4(0,blocksize)

      Pspb.new(index + 1,entry.definition,pspb1,pspb2,pspb3,pspb4)
    end

    puts "node1\tnode2"
    pspbs.combination(2).each do |i|
      puts "#{i[0].name}\t#{i[1].name}" if (i[0].all_pspbs & i[1].all_pspbs).size > 0
    end
  end

  command :attribute do |a|
    a.syntax = 'network attribute [options]'
    a.summary = 'generate network attributes file'
    a.description = 'This command generates a network attributes file for loading into cystoscape'

    a.example 'network', 'attribute --infile FILENAME'
    a.option '-f','--infile FILE', 'specify the sequence input file'

    #a structure to hold the cyspolv classification and bs sharing groups
    class Attr < Struct.new(:id,:name,:cyspolv,:bsgroup);end

    a.action do |args, options|
      seqfile = options.infile

      attrs = Bio::FlatFile.open(seqfile).each_with_index.map do |entry,index|
        cyspolv = Bio::Sequence::AA.new(entry.seq).cyspolv_group
        bsgroup = Bio::Sequence::AA.new(entry.seq).bs_group
        Attr.new(index + 1,entry.definition,cyspolv,bsgroup)
      end

      puts "gene\tcyspolvgroup\tbsgroup"
      attrs.each do |attri|
        puts "#{attri.name}\t#{attri.cyspolv}\t#{attri.bsgroup}"
      end
    end
  end
end
