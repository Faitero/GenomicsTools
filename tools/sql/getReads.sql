select distinct Sample.number, organism, concat(propertyValue, year(Request.createDate), '/', filename), codeApplication, codeSampleFileType
    from Request join Sample on (Sample.idRequest = Request.idRequest)
    join Organism on Organism.idOrganism = Sample.idOrganism
    join SampleExperimentFile on SampleExperimentFile.idSample = Sample.idSample
    join ExperimentFile on ExperimentFile.idExperimentFile = SampleExperimentFile.idExperimentFile
    join PropertyDictionary on PropertyDictionary.propertyName='experiment_directory'
    where SampleExperimentFile.codeSampleFileType = 'fastqRead1' or SampleExperimentFile.codeSampleFileType = 'fastqRead2';
