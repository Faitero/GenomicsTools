#!/usr/bin/perl -w

use strict;
use b2b::db;

my $dbh = b2b::db::dbConnect();

b2b::db::buildDatabase( dbh=>$dbh );