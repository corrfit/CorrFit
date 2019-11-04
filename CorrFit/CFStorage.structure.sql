-- MySQL dump 9.09
--
-- Host: localhost    Database: CFStorage
---------------------------------------------------------
-- Server version	4.0.15

--
-- Table structure for table `FitParameters`
--

CREATE TABLE FitParameters (
  fitParameterID int(11) NOT NULL auto_increment,
  fitBinsCount smallint(6) NOT NULL default '0',
  fitMin varchar(9) NOT NULL default '0.000',
  fitMax varchar(9) NOT NULL default '0.000',
  normMin varchar(9) NOT NULL default '0.000',
  normMax varchar(9) NOT NULL default '0.000',
  momResScale varchar(9) NOT NULL default '0.000',
  PRIMARY KEY  (fitParameterID),
  KEY fitBinsCount (fitBinsCount,fitMin,fitMax,normMin,normMax,momResScale)
) TYPE=MyISAM;

--
-- Table structure for table `FunctionData`
--

CREATE TABLE FunctionData (
  functionDataID bigint(20) NOT NULL auto_increment,
  modelParameterID bigint(20) NOT NULL default '0',
  pairHashID int(11) NOT NULL default '0',
  fitParameterID int(11) NOT NULL default '0',
  interactionID smallint(6) NOT NULL default '0',
  pairTypeID smallint(6) NOT NULL default '0',
  data text NOT NULL,
  PRIMARY KEY  (functionDataID),
  UNIQUE KEY FD_FunctionDataUniqe (pairTypeID,interactionID,fitParameterID,pairHashID,modelParameterID,data(5)),
  KEY FD_ModelParametersIndex (modelParameterID),
  KEY FD_checkData (modelParameterID,pairHashID,fitParameterID,interactionID,pairTypeID)
) TYPE=MyISAM;

--
-- Table structure for table `Interactions`
--

CREATE TABLE Interactions (
  interactionID smallint(6) NOT NULL auto_increment,
  strong tinyint(1) NOT NULL default '0',
  coulomb tinyint(1) NOT NULL default '0',
  quantumStatistics tinyint(1) NOT NULL default '0',
  PRIMARY KEY  (interactionID),
  UNIQUE KEY IN_InteractionUnique (strong,coulomb,quantumStatistics)
) TYPE=MyISAM;

--
-- Table structure for table `ModelParameters`
--

CREATE TABLE ModelParameters (
  modelParameterID bigint(20) NOT NULL auto_increment,
  sourceModelID smallint(6) NOT NULL default '0',
  parameter1 varchar(9) NOT NULL default '0.000',
  parameter2 varchar(9) NOT NULL default '0.000',
  parameter3 varchar(9) NOT NULL default '0.000',
  parameter4 varchar(9) NOT NULL default '0.000',
  parameter5 varchar(9) NOT NULL default '0.000',
  parameter6 varchar(9) NOT NULL default '0.000',
  parameter7 varchar(9) NOT NULL default '0.000',
  parameter8 varchar(9) NOT NULL default '0.000',
  parameter9 varchar(9) NOT NULL default '0.000',
  parameter10 varchar(9) NOT NULL default '0.000',
  PRIMARY KEY  (modelParameterID),
  UNIQUE KEY MP_ModelParametersUnique (sourceModelID,parameter1,parameter2,parameter3,parameter4,parameter5,parameter6,parameter7,parameter8,parameter9,parameter10)
) TYPE=MyISAM;

--
-- Table structure for table `PairHashes`
--

CREATE TABLE PairHashes (
  pairHashID int(11) NOT NULL auto_increment,
  sampleHash varchar(30) NOT NULL default '',
  pairCount int(11) NOT NULL default '0',
  PRIMARY KEY  (pairHashID),
  UNIQUE KEY PH_PairHashUnique (sampleHash,pairCount)
) TYPE=MyISAM;

--
-- Table structure for table `PairTypes`
--

CREATE TABLE PairTypes (
  pairTypeID smallint(6) NOT NULL auto_increment,
  pairTypeName varchar(30) NOT NULL default '',
  pairTypeSymbol varchar(30) NOT NULL default '',
  PRIMARY KEY  (pairTypeID),
  UNIQUE KEY pairTypeName (pairTypeName)
) TYPE=MyISAM;

--
-- Table structure for table `SourceModels`
--

CREATE TABLE SourceModels (
  sourceModelID smallint(6) NOT NULL auto_increment,
  sourceModelName varchar(30) NOT NULL default '',
  parameterCount smallint(6) NOT NULL default '0',
  PRIMARY KEY  (sourceModelID),
  UNIQUE KEY SM_SourceModelUnique (sourceModelName,parameterCount)
) TYPE=MyISAM;

