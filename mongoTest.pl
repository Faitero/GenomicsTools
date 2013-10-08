use strict;
use warnings;
use MongoDB;
use MongoDB::OID;
use Data::Dumper;
use Tie::IxHash;
use tools::meta;
use tools::mongo;


my 	 $samplesheet= "/home/laurentt/sampleSheet/gnomex_sampleannotations_SEIDMAN_030713.txt";
my 	 $alignGTF;						
my	 $dataPath = "/Data01/gnomex/ExperimentData/";
my	 $analysisDir = "/Data01/gnomex/Analysis/experiment/";
my $exp = "166R";

my $sampleHash = tools::meta::parseSampleSheet($samplesheet);

# the following tests the mongo.pm addSamps sub

tools::mongo::addSamps(sh=> $sampleHash);

tools::mongo::addBamStat(exp => "77R");

my $bamstatHash = tools::mongo::getBamStat(exp => "77R");

# print "bamstatHash".Dumper($bamstatHash)."\n";

# for my $samp (keys %$bamstatHash){
#     print "samp\t$samp\n";
#     for my $fld ( keys %{$bamstatHash->{$samp}->{data}} ){
#         print "$fld\t".$bamstatHash->{$samp}->{data}->{$fld}."\n";
#     }
# }

tools::mongo::addMarkDuplicates(exp => "77R");

tools::mongo::addInnerDistance(exp => "77R");

tools::mongo::addGeneBodyCov(exp => "77R");
#  uncomment to do the experiment on only one sample
# my $sampleHash = tools::meta::makeExpSampleHash(
# 	exp => $exp,
# 	samplesheet => $samplesheet,
# 	);



# my $client = MongoDB::MongoClient->new;

# my $db = $client->get_database( 'b2btest' );

# my $experiment = $db->get_collection( 'experiment' );

#  $experiment->drop;

# for my $key (keys %$sampleHash){
#     my $expAtts = $sampleHash->{$key};
#     print "key\t$key\n";
#     print Dumper($expAtts);
#     $experiment->update({"_id"=>$key}, { '$set' => {%$expAtts}},{"upsert"=> 1} );
# }
 

#  my $expColl = $db->get_collection('exptest');

#  print "exp\t$exp\n";

#  $expColl->update({"_id"=>$exp}, {'$set' => { } },{"upsert"=> 1});

# #my $id = $users->insert({"name"=>"Bill"});

# #print "id\t$id\n";

# #$users->remove({"name" => "Joe"});

# my $experiments = $experiment->find->fields({});

# while ( my $exp = $experiments->next){
#     print $exp->{'_id'}."\n";
# }




# print Dumper($experiment)."\n";

#my $some_users = $users->find({"name"=>"Joe"});

#print $some_users->{'name'}."\n";

#print Dumper($some_users)."\n";


#while (my $doc = $experiments->next){
    #print $doc->{'id'}."\n";
   #}

   ## to use run command use the Tie::IxHash  object --- it is an ordered Hashi

#my $cmd = Tie::IxHash->new("create" => "posts""
        #"capped" => boolean::true,
            #"size" => 10240,
                #"max" => 100);

                #$db->run_command($cmd);



