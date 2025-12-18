use std::path::PathBuf;
use std::str::FromStr;

use ahash::HashMap;
use anyhow::Result;
use itertools::izip;
use ndarray::Array2;
use noodles::core::{Position, Region};
use noodles::cram::record::feature;
use noodles::{bam, bed};
use polars::prelude::*;
use rayon::prelude::*;

use crate::bam_utils::{BamStats, Iv, regions_to_lapper};

// Aim: to convert the following python code to rust

// class heatmapper(object):
//     """
//     Class to handle the reading and
//     plotting of matrices.
//     """

//     def __init__(self):
//         self.parameters = None
//         self.lengthDict = None
//         self.matrix = None
//         self.regions = None
//         self.blackList = None
//         self.quiet = True
//         # These are parameters that were single values in versions <3 but are now internally lists. See issue #614
//         self.special_params = set(['unscaled 5 prime', 'unscaled 3 prime', 'body', 'downstream', 'upstream', 'ref point', 'bin size'])

//     def getTicks(self, idx):
//         """
//         This is essentially a wrapper around getProfileTicks to accomdate the fact that each column has its own ticks.
//         """
//         xticks, xtickslabel = getProfileTicks(self, self.reference_point_label[idx], self.startLabel, self.endLabel, idx)
//         return xticks, xtickslabel

//     def computeMatrix(self, score_file_list, regions_file, parameters, blackListFileName=None, verbose=False, allArgs=None):
//         """
//         Splits into
//         multiple cores the computation of the scores
//         per bin for each region (defined by a hash '#'
//         in the regions (BED/GFF) file.
//         """
//         if parameters['body'] > 0 and \
//                 parameters['body'] % parameters['bin size'] > 0:
//             exit("The --regionBodyLength has to be "
//                  "a multiple of --binSize.\nCurrently the "
//                  "values are {} {} for\nregionsBodyLength and "
//                  "binSize respectively\n".format(parameters['body'],
//                                                  parameters['bin size']))

//         # the beforeRegionStartLength is extended such that
//         # length is a multiple of binSize
//         if parameters['downstream'] % parameters['bin size'] > 0:
//             exit("Length of region after the body has to be "
//                  "a multiple of --binSize.\nCurrent value "
//                  "is {}\n".format(parameters['downstream']))

//         if parameters['upstream'] % parameters['bin size'] > 0:
//             exit("Length of region before the body has to be a multiple of "
//                  "--binSize\nCurrent value is {}\n".format(parameters['upstream']))

//         if parameters['unscaled 5 prime'] % parameters['bin size'] > 0:
//             exit("Length of the unscaled 5 prime region has to be a multiple of "
//                  "--binSize\nCurrent value is {}\n".format(parameters['unscaled 5 prime']))

//         if parameters['unscaled 3 prime'] % parameters['bin size'] > 0:
//             exit("Length of the unscaled 5 prime region has to be a multiple of "
//                  "--binSize\nCurrent value is {}\n".format(parameters['unscaled 3 prime']))

//         if parameters['unscaled 5 prime'] + parameters['unscaled 3 prime'] > 0 and parameters['body'] == 0:
//             exit('Unscaled 5- and 3-prime regions only make sense with the scale-regions subcommand.\n')

//         # Take care of GTF options
//         transcriptID = "transcript"
//         exonID = "exon"
//         transcript_id_designator = "transcript_id"
//         keepExons = False
//         self.quiet = False
//         if allArgs is not None:
//             allArgs = vars(allArgs)
//             transcriptID = allArgs.get("transcriptID", transcriptID)
//             exonID = allArgs.get("exonID", exonID)
//             transcript_id_designator = allArgs.get("transcript_id_designator", transcript_id_designator)
//             keepExons = allArgs.get("keepExons", keepExons)
//             self.quiet = allArgs.get("quiet", self.quiet)

//         chromSizes, _ = getScorePerBigWigBin.getChromSizes(score_file_list)
//         res, labels = mapReduce.mapReduce([score_file_list, parameters],
//                                           compute_sub_matrix_wrapper,
//                                           chromSizes,
//                                           self_=self,
//                                           bedFile=regions_file,
//                                           blackListFileName=blackListFileName,
//                                           numberOfProcessors=parameters['proc number'],
//                                           includeLabels=True,
//                                           transcriptID=transcriptID,
//                                           exonID=exonID,
//                                           transcript_id_designator=transcript_id_designator,
//                                           keepExons=keepExons,
//                                           verbose=verbose)
//         # each worker in the pool returns a tuple containing
//         # the submatrix data, the regions that correspond to the
//         # submatrix, and the number of regions lacking scores
//         # Since this is largely unsorted, we need to sort by group

//         # merge all the submatrices into matrix
//         matrix = np.concatenate([r[0] for r in res], axis=0)
//         regions = []
//         regions_no_score = 0
//         for idx in range(len(res)):
//             if len(res[idx][1]):
//                 regions.extend(res[idx][1])
//                 regions_no_score += res[idx][2]
//         groups = [x[3] for x in regions]
//         foo = sorted(zip(groups, list(range(len(regions))), regions))
//         sortIdx = [x[1] for x in foo]
//         regions = [x[2] for x in foo]
//         matrix = matrix[sortIdx]

//         # mask invalid (nan) values
//         matrix = np.ma.masked_invalid(matrix)

//         assert matrix.shape[0] == len(regions), \
//             "matrix length does not match regions length"

//         if len(regions) == 0:
//             sys.stderr.write("\nERROR: Either the BED file does not contain any valid regions or there are none remaining after filtering.\n")
//             exit(1)
//         if regions_no_score == len(regions):
//             exit("\nERROR: None of the BED regions could be found in the bigWig"
//                  "file.\nPlease check that the bigwig file is valid and "
//                  "that the chromosome names between the BED file and "
//                  "the bigWig file correspond to each other\n")

//         if regions_no_score > len(regions) * 0.75:
//             file_type = 'bigwig' if score_file_list[0].endswith(".bw") else "BAM"
//             prcnt = 100 * float(regions_no_score) / len(regions)
//             sys.stderr.write(
//                 "\n\nWarning: {0:.2f}% of regions are *not* associated\n"
//                 "to any score in the given {1} file. Check that the\n"
//                 "chromosome names from the BED file are consistent with\n"
//                 "the chromosome names in the given {2} file and that both\n"
//                 "files refer to the same species\n\n".format(prcnt,
//                                                              file_type,
//                                                              file_type))

//         self.parameters = parameters

//         numcols = matrix.shape[1]
//         num_ind_cols = self.get_num_individual_matrix_cols()
//         sample_boundaries = list(range(0, numcols + num_ind_cols, num_ind_cols))
//         if allArgs is not None and allArgs['samplesLabel'] is not None:
//             sample_labels = allArgs['samplesLabel']
//         else:
//             sample_labels = smartLabels(score_file_list)

//         # Determine the group boundaries
//         group_boundaries = []
//         group_labels_filtered = []
//         last_idx = -1
//         for x in range(len(regions)):
//             if regions[x][3] != last_idx:
//                 last_idx = regions[x][3]
//                 group_boundaries.append(x)
//                 group_labels_filtered.append(labels[last_idx])
//         group_boundaries.append(len(regions))

//         # check if a given group is too small. Groups that
//         # are too small can't be plotted and an exception is thrown.
//         group_len = np.diff(group_boundaries)
//         if len(group_len) > 1:
//             sum_len = sum(group_len)
//             group_frac = [float(x) / sum_len for x in group_len]
//             if min(group_frac) <= 0.002:
//                 sys.stderr.write(
//                     "One of the groups defined in the bed file is "
//                     "too small.\nGroups that are too small can't be plotted. "
//                     "\n")

//         self.matrix = _matrix(regions, matrix,
//                               group_boundaries,
//                               sample_boundaries,
//                               group_labels_filtered,
//                               sample_labels)

//         if parameters['skip zeros']:
//             self.matrix.removeempty()

//     @staticmethod
//     def compute_sub_matrix_worker(self, chrom, start, end, score_file_list, parameters, regions):
//         """
//         Returns
//         -------
//         numpy matrix
//             A numpy matrix that contains per each row the values found per each of the regions given
//         """
//         if parameters['verbose']:
//             sys.stderr.write("Processing {}:{}-{}\n".format(chrom, start, end))

//         # read BAM or scores file
//         score_file_handles = []
//         for sc_file in score_file_list:
//             score_file_handles.append(pyBigWig.open(sc_file))

//         # determine the number of matrix columns based on the lengths
//         # given by the user, times the number of score files
//         matrix_cols = len(score_file_list) * \
//             ((parameters['downstream'] +
//               parameters['unscaled 5 prime'] + parameters['unscaled 3 prime'] +
//               parameters['upstream'] + parameters['body']) //
//              parameters['bin size'])

//         # create an empty matrix to store the values
//         sub_matrix = np.zeros((len(regions), matrix_cols))
//         sub_matrix[:] = np.NAN

//         j = 0
//         sub_regions = []
//         regions_no_score = 0
//         for transcript in regions:
//             feature_chrom = transcript[0]
//             exons = transcript[1]
//             feature_start = exons[0][0]
//             feature_end = exons[-1][1]
//             feature_name = transcript[2]
//             feature_strand = transcript[4]
//             padLeft = 0
//             padRight = 0
//             padLeftNaN = 0
//             padRightNaN = 0
//             upstream = []
//             downstream = []

//             # get the body length
//             body_length = np.sum([x[1] - x[0] for x in exons]) - parameters['unscaled 5 prime'] - parameters['unscaled 3 prime']

//             # print some information
//             if parameters['body'] > 0 and \
//                     body_length < parameters['bin size']:
//                 if not self.quiet:
//                     sys.stderr.write("A region that is shorter than the bin size (possibly only after accounting for unscaled regions) was found: "
//                                      "({0}) {1} {2}:{3}:{4}. Skipping...\n".format((body_length - parameters['unscaled 5 prime'] - parameters['unscaled 3 prime']),
//                                                                                    feature_name, feature_chrom,
//                                                                                    feature_start, feature_end))
//                 coverage = np.zeros(matrix_cols)
//                 if not parameters['missing data as zero']:
//                     coverage[:] = np.nan
//             else:
//                 if feature_strand == '-':
//                     if parameters['downstream'] > 0:
//                         upstream = [(feature_start - parameters['downstream'], feature_start)]
//                     if parameters['upstream'] > 0:
//                         downstream = [(feature_end, feature_end + parameters['upstream'])]
//                     unscaled5prime, body, unscaled3prime, padLeft, padRight = chopRegions(exons, left=parameters['unscaled 3 prime'], right=parameters['unscaled 5 prime'])
//                     # bins per zone
//                     a = parameters['downstream'] // parameters['bin size']
//                     b = parameters['unscaled 3 prime'] // parameters['bin size']
//                     d = parameters['unscaled 5 prime'] // parameters['bin size']
//                     e = parameters['upstream'] // parameters['bin size']
//                 else:
//                     if parameters['upstream'] > 0:
//                         upstream = [(feature_start - parameters['upstream'], feature_start)]
//                     if parameters['downstream'] > 0:
//                         downstream = [(feature_end, feature_end + parameters['downstream'])]
//                     unscaled5prime, body, unscaled3prime, padLeft, padRight = chopRegions(exons, left=parameters['unscaled 5 prime'], right=parameters['unscaled 3 prime'])
//                     a = parameters['upstream'] // parameters['bin size']
//                     b = parameters['unscaled 5 prime'] // parameters['bin size']
//                     d = parameters['unscaled 3 prime'] // parameters['bin size']
//                     e = parameters['downstream'] // parameters['bin size']
//                 c = parameters['body'] // parameters['bin size']

//                 # build zones (each is a list of tuples)
//                 #  zone0: region before the region start,
//                 #  zone1: unscaled 5 prime region
//                 #  zone2: the body of the region
//                 #  zone3: unscaled 3 prime region
//                 #  zone4: the region from the end of the region downstream
//                 #  the format for each zone is: [(start, end), ...], number of bins
//                 # Note that for "reference-point", upstream/downstream will go
//                 # through the exons (if requested) and then possibly continue
//                 # on the other side (unless parameters['nan after end'] is true)
//                 if parameters['body'] > 0:
//                     zones = [(upstream, a), (unscaled5prime, b), (body, c), (unscaled3prime, d), (downstream, e)]
//                 elif parameters['ref point'] == 'TES':  # around TES
//                     if feature_strand == '-':
//                         downstream, body, unscaled3prime, padRight, _ = chopRegions(exons, left=parameters['upstream'])
//                         if padRight > 0 and parameters['nan after end'] is True:
//                             padRightNaN += padRight
//                         elif padRight > 0:
//                             downstream.append((downstream[-1][1], downstream[-1][1] + padRight))
//                         padRight = 0
//                     else:
//                         unscale5prime, body, upstream, _, padLeft = chopRegions(exons, right=parameters['upstream'])
//                         if padLeft > 0 and parameters['nan after end'] is True:
//                             padLeftNaN += padLeft
//                         elif padLeft > 0:
//                             upstream.insert(0, (upstream[0][0] - padLeft, upstream[0][0]))
//                         padLeft = 0
//                     e = np.sum([x[1] - x[0] for x in downstream]) // parameters['bin size']
//                     a = np.sum([x[1] - x[0] for x in upstream]) // parameters['bin size']
//                     zones = [(upstream, a), (downstream, e)]
//                 elif parameters['ref point'] == 'center':  # at the region center
//                     if feature_strand == '-':
//                         upstream, downstream, padLeft, padRight = chopRegionsFromMiddle(exons, left=parameters['downstream'], right=parameters['upstream'])
//                     else:
//                         upstream, downstream, padLeft, padRight = chopRegionsFromMiddle(exons, left=parameters['upstream'], right=parameters['downstream'])
//                     if padLeft > 0 and parameters['nan after end'] is True:
//                         padLeftNaN += padLeft
//                     elif padLeft > 0:
//                         if len(upstream) > 0:
//                             upstream.insert(0, (upstream[0][0] - padLeft, upstream[0][0]))
//                         else:
//                             upstream = [(downstream[0][0] - padLeft, downstream[0][0])]
//                     padLeft = 0
//                     if padRight > 0 and parameters['nan after end'] is True:
//                         padRightNaN += padRight
//                     elif padRight > 0:
//                         downstream.append((downstream[-1][1], downstream[-1][1] + padRight))
//                     padRight = 0
//                     a = np.sum([x[1] - x[0] for x in upstream]) // parameters['bin size']
//                     e = np.sum([x[1] - x[0] for x in downstream]) // parameters['bin size']
//                     # It's possible for a/e to be floats or 0 yet upstream/downstream isn't empty
//                     if a < 1:
//                         upstream = []
//                         a = 0
//                     if e < 1:
//                         downstream = []
//                         e = 0
//                     zones = [(upstream, a), (downstream, e)]
//                 else:  # around TSS
//                     if feature_strand == '-':
//                         unscale5prime, body, upstream, _, padLeft = chopRegions(exons, right=parameters['downstream'])
//                         if padLeft > 0 and parameters['nan after end'] is True:
//                             padLeftNaN += padLeft
//                         elif padLeft > 0:
//                             upstream.insert(0, (upstream[0][0] - padLeft, upstream[0][0]))
//                         padLeft = 0
//                     else:
//                         downstream, body, unscaled3prime, padRight, _ = chopRegions(exons, left=parameters['downstream'])
//                         if padRight > 0 and parameters['nan after end'] is True:
//                             padRightNaN += padRight
//                         elif padRight > 0:
//                             downstream.append((downstream[-1][1], downstream[-1][1] + padRight))
//                         padRight = 0
//                     a = np.sum([x[1] - x[0] for x in upstream]) // parameters['bin size']
//                     e = np.sum([x[1] - x[0] for x in downstream]) // parameters['bin size']
//                     zones = [(upstream, a), (downstream, e)]

//                 foo = parameters['upstream']
//                 bar = parameters['downstream']
//                 if feature_strand == '-':
//                     foo, bar = bar, foo
//                 if padLeftNaN > 0:
//                     expected = foo // parameters['bin size']
//                     padLeftNaN = int(round(float(padLeftNaN) / parameters['bin size']))
//                     if expected - padLeftNaN - a > 0:
//                         padLeftNaN += 1
//                 if padRightNaN > 0:
//                     expected = bar // parameters['bin size']
//                     padRightNaN = int(round(float(padRightNaN) / parameters['bin size']))
//                     if expected - padRightNaN - e > 0:
//                         padRightNaN += 1

//                 coverage = []
//                 # compute the values for each of the files being processed.
//                 # "cov" is a numpy array of bins
//                 for sc_handler in score_file_handles:
//                     # We're only supporting bigWig files at this point
//                     cov = heatmapper.coverage_from_big_wig(
//                         sc_handler, feature_chrom, zones,
//                         parameters['bin size'],
//                         parameters['bin avg type'],
//                         parameters['missing data as zero'],
//                         not self.quiet)

//                     if padLeftNaN > 0:
//                         cov = np.concatenate([[np.nan] * padLeftNaN, cov])
//                     if padRightNaN > 0:
//                         cov = np.concatenate([cov, [np.nan] * padRightNaN])

//                     if feature_strand == "-":
//                         cov = cov[::-1]

//                     coverage = np.hstack([coverage, cov])

//             if coverage is None:
//                 regions_no_score += 1
//                 if not self.quiet:
//                     sys.stderr.write(
//                         "No data was found for region "
//                         "{0} {1}:{2}-{3}. Skipping...\n".format(
//                             feature_name, feature_chrom,
//                             feature_start, feature_end))

//                 coverage = np.zeros(matrix_cols)
//                 if not parameters['missing data as zero']:
//                     coverage[:] = np.nan

//             try:
//                 temp = coverage.copy()
//                 temp[np.isnan(temp)] = 0
//             except:
//                 if not self.quiet:
//                     sys.stderr.write(
//                         "No scores defined for region "
//                         "{0} {1}:{2}-{3}. Skipping...\n".format(feature_name,
//                                                                 feature_chrom,
//                                                                 feature_start,
//                                                                 feature_end))
//                 coverage = np.zeros(matrix_cols)
//                 if not parameters['missing data as zero']:
//                     coverage[:] = np.nan

//             if parameters['min threshold'] is not None and coverage.min() <= parameters['min threshold']:
//                 continue
//             if parameters['max threshold'] is not None and coverage.max() >= parameters['max threshold']:
//                 continue
//             if parameters['scale'] != 1:
//                 coverage = parameters['scale'] * coverage

//             sub_matrix[j, :] = coverage

//             sub_regions.append(transcript)
//             j += 1

//         # remove empty rows
//         sub_matrix = sub_matrix[0:j, :]
//         if len(sub_regions) != len(sub_matrix[:, 0]):
//             sys.stderr.write("regions lengths do not match\n")
//         return sub_matrix, sub_regions, regions_no_score

//     @staticmethod
//     def coverage_from_array(valuesArray, zones, binSize, avgType):
//         try:
//             valuesArray[0]
//         except (IndexError, TypeError) as detail:
//             sys.stderr.write("{0}\nvalues array value: {1}, zones {2}\n".format(detail, valuesArray, zones))

//         cvglist = []
//         zoneEnd = 0
//         valStart = 0
//         valEnd = 0
//         for zone, nBins in zones:
//             if nBins:
//                 # linspace is used to more or less evenly partition the data points into the given number of bins
//                 zoneEnd += nBins
//                 valStart = valEnd
//                 valEnd += np.sum([x[1] - x[0] for x in zone])
//                 counts_list = []

//                 # Partition the space into bins
//                 if nBins == 1:
//                     pos_array = np.array([valStart])
//                 else:
//                     pos_array = np.linspace(valStart, valEnd, nBins, endpoint=False, dtype=int)
//                 pos_array = np.append(pos_array, valEnd)

//                 idx = 0
//                 while idx < nBins:
//                     idxStart = int(pos_array[idx])
//                     idxEnd = max(int(pos_array[idx + 1]), idxStart + 1)
//                     try:
//                         counts_list.append(heatmapper.my_average(valuesArray[idxStart:idxEnd], avgType))
//                     except Exception as detail:
//                         sys.stderr.write("Exception found: {0}\n".format(detail))
//                     idx += 1
//                 cvglist.append(np.array(counts_list))

//         return np.concatenate(cvglist)

//     @staticmethod
//     def change_chrom_names(chrom):
//         """
//         Changes UCSC chromosome names to ensembl chromosome names
//         and vice versa.
//         """
//         if chrom.startswith('chr'):
//             # remove the chr part from chromosome name
//             chrom = chrom[3:]
//             if chrom == "M":
//                 chrom = "MT"
//         else:
//             # prefix with 'chr' the chromosome name
//             chrom = 'chr' + chrom
//             if chrom == "chrMT":
//                 chrom = "chrM"

//         return chrom

//     @staticmethod
//     def coverage_from_big_wig(bigwig, chrom, zones, binSize, avgType, nansAsZeros=False, verbose=True):

//         """
//         uses pyBigWig
//         to query a region define by chrom and zones.
//         The output is an array that contains the bigwig
//         value per base pair. The summary over bins is
//         done in a later step when coverage_from_array is called.
//         This method is more reliable than querying the bins
//         directly from the bigwig, which should be more efficient.

//         By default, any region, even if no chromosome match is found
//         on the bigwig file, produces a result. In other words
//         no regions are skipped.

//         zones: array as follows zone0: region before the region start,
//                                 zone1: 5' unscaled region (if present)
//                                 zone2: the body of the region (not always present)
//                                 zone3: 3' unscaled region (if present)
//                                 zone4: the region from the end of the region downstream

//                each zone is a tuple containing start, end, and number of bins

//         This is useful if several matrices wants to be merged
//         or if the sorted BED output of one computeMatrix operation
//         needs to be used for other cases
//         """
//         nVals = 0
//         for zone, _ in zones:
//             for region in zone:
//                 nVals += region[1] - region[0]

//         values_array = np.zeros(nVals)
//         if not nansAsZeros:
//             values_array[:] = np.nan
//         if chrom not in list(bigwig.chroms().keys()):
//             unmod_name = chrom
//             chrom = heatmapper.change_chrom_names(chrom)
//             if chrom not in list(bigwig.chroms().keys()):
//                 if verbose:
//                     sys.stderr.write("Warning: Your chromosome names do not match.\nPlease check that the "
//                                      "chromosome names in your BED file\ncorrespond to the names in your "
//                                      "bigWig file.\nAn empty line will be added to your heatmap.\nThe problematic "
//                                      "chromosome name is {0}\n\n".format(unmod_name))

//                 # return empty nan array
//                 return heatmapper.coverage_from_array(values_array, zones, binSize, avgType)

//         maxLen = bigwig.chroms(chrom)
//         startIdx = 0
//         endIdx = 0
//         for zone, _ in zones:
//             for region in zone:
//                 startIdx = endIdx
//                 if region[0] < 0:
//                     endIdx += abs(region[0])
//                     values_array[startIdx:endIdx] = np.nan
//                     startIdx = endIdx
//                 start = max(0, region[0])
//                 end = min(maxLen, region[1])
//                 endIdx += end - start
//                 if start < end:
//                     # This won't be the case if we extend off the front of a chromosome, such as (-100, 0)
//                     values_array[startIdx:endIdx] = bigwig.values(chrom, start, end)
//                 if end < region[1]:
//                     startIdx = endIdx
//                     endIdx += region[1] - end
//                     values_array[startIdx:endIdx] = np.nan

//         # replaces nans for zeros
//         if nansAsZeros:
//             values_array[np.isnan(values_array)] = 0

//         return heatmapper.coverage_from_array(values_array, zones,
//                                               binSize, avgType)

//     @staticmethod
//     def my_average(valuesArray, avgType='mean'):
//         """
//         computes the mean, median, etc but only for those values
//         that are not Nan
//         """
//         valuesArray = np.ma.masked_invalid(valuesArray)
//         avg = np.ma.__getattribute__(avgType)(valuesArray)
//         if isinstance(avg, np.ma.core.MaskedConstant):
//             return np.nan
//         else:
//             return avg

//     def matrix_from_dict(self, matrixDict, regionsDict, parameters):
//         self.regionsDict = regionsDict
//         self.matrixDict = matrixDict
//         self.parameters = parameters
//         self.lengthDict = OrderedDict()
//         self.matrixAvgsDict = OrderedDict()

//     def read_matrix_file(self, matrix_file):
//         # reads a bed file containing the position
//         # of genomic intervals
//         # In case a hash sign '#' is found in the
//         # file, this is considered as a delimiter
//         # to split the heatmap into groups

//         import json
//         regions = []
//         matrix_rows = []
//         current_group_index = 0
//         max_group_bound = None

//         fh = gzip.open(matrix_file)
//         for line in fh:
//             line = toString(line).strip()
//             # read the header file containing the parameters
//             # used
//             if line.startswith("@"):
//                 # the parameters used are saved using
//                 # json
//                 self.parameters = json.loads(line[1:].strip())
//                 max_group_bound = self.parameters['group_boundaries'][1]
//                 continue

//             # split the line into bed interval and matrix values
//             region = line.split('\t')
//             chrom, start, end, name, score, strand = region[0:6]
//             matrix_row = np.ma.masked_invalid(np.fromiter(region[6:], float))
//             matrix_rows.append(matrix_row)
//             starts = start.split(",")
//             ends = end.split(",")
//             regs = [(int(x), int(y)) for x, y in zip(starts, ends)]
//             # get the group index
//             if len(regions) >= max_group_bound:
//                 current_group_index += 1
//                 max_group_bound = self.parameters['group_boundaries'][current_group_index + 1]
//             regions.append([chrom, regs, name, max_group_bound, strand, score])

//         matrix = np.vstack(matrix_rows)
//         self.matrix = _matrix(regions, matrix, self.parameters['group_boundaries'],
//                               self.parameters['sample_boundaries'],
//                               group_labels=self.parameters['group_labels'],
//                               sample_labels=self.parameters['sample_labels'])

//         if 'sort regions' in self.parameters:
//             self.matrix.set_sorting_method(self.parameters['sort regions'],
//                                            self.parameters['sort using'])

//         # Versions of computeMatrix before 3.0 didn't have an entry of these per column, fix that
//         nSamples = len(self.matrix.sample_labels)
//         h = dict()
//         for k, v in self.parameters.items():
//             if k in self.special_params and type(v) is not list:
//                 v = [v] * nSamples
//                 if len(v) == 0:
//                     v = [None] * nSamples
//             h[k] = v
//         self.parameters = h

//         return

//     def save_matrix(self, file_name):
//         """
//         saves the data required to reconstruct the matrix
//         the format is:
//         A header containing the parameters used to create the matrix
//         encoded as:
//         @key:value\tkey2:value2 etc...
//         The rest of the file has the same first 5 columns of a
//         BED file: chromosome name, start, end, name, score and strand,
//         all separated by tabs. After the fifth column the matrix
//         values are appended separated by tabs.
//         Groups are separated by adding a line starting with a hash (#)
//         and followed by the group name.

//         The file is gzipped.
//         """
//         import json
//         self.parameters['sample_labels'] = self.matrix.sample_labels
//         self.parameters['group_labels'] = self.matrix.group_labels
//         self.parameters['sample_boundaries'] = self.matrix.sample_boundaries
//         self.parameters['group_boundaries'] = self.matrix.group_boundaries

//         # Redo the parameters, ensuring things related to ticks and labels are repeated appropriately
//         nSamples = len(self.matrix.sample_labels)
//         h = dict()
//         for k, v in self.parameters.items():
//             if type(v) is list and len(v) == 0:
//                 v = None
//             if k in self.special_params and type(v) is not list:
//                 v = [v] * nSamples
//                 if len(v) == 0:
//                     v = [None] * nSamples
//             h[k] = v
//         fh = gzip.open(file_name, 'wb')
//         params_str = json.dumps(h, separators=(',', ':'))
//         fh.write(toBytes("@" + params_str + "\n"))
//         score_list = np.ma.masked_invalid(np.mean(self.matrix.matrix, axis=1))
//         for idx, region in enumerate(self.matrix.regions):
//             # join np_array values
//             # keeping nans while converting them to strings
//             if not np.ma.is_masked(score_list[idx]):
//                 float(score_list[idx])
//             matrix_values = "\t".join(
//                 np.char.mod('%f', self.matrix.matrix[idx, :]))
//             starts = ["{0}".format(x[0]) for x in region[1]]
//             ends = ["{0}".format(x[1]) for x in region[1]]
//             starts = ",".join(starts)
//             ends = ",".join(ends)
//             # BEDish format (we don't currently store the score)
//             fh.write(
//                 toBytes('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(
//                         region[0],
//                         starts,
//                         ends,
//                         region[2],
//                         region[5],
//                         region[4],
//                         matrix_values)))
//         fh.close()

//     def save_tabulated_values(self, file_handle, reference_point_label='TSS', start_label='TSS', end_label='TES', averagetype='mean'):
//         """
//         Saves the values averaged by col using the avg_type
//         given

//         Args:
//             file_handle: file name to save the file
//             reference_point_label: Name of the reference point label
//             start_label: Name of the star label
//             end_label: Name of the end label
//             averagetype: average type (e.g. mean, median, std)

//         """
//         #  get X labels
//         w = self.parameters['bin size']
//         b = self.parameters['upstream']
//         a = self.parameters['downstream']
//         c = self.parameters.get('unscaled 5 prime', 0)
//         d = self.parameters.get('unscaled 3 prime', 0)
//         m = self.parameters['body']

//         xticks = []
//         xtickslabel = []
//         for idx in range(self.matrix.get_num_samples()):
//             if b[idx] < 1e5:
//                 quotient = 1000
//                 symbol = 'Kb'
//             else:
//                 quotient = 1e6
//                 symbol = 'Mb'

//             if m[idx] == 0:
//                 last = 0
//                 if len(xticks):
//                     last = xticks[-1]
//                 xticks.extend([last + (k / w[idx]) for k in [w[idx], b[idx], b[idx] + a[idx]]])
//                 xtickslabel.extend(['{0:.1f}{1}'.format(-(float(b[idx]) / quotient), symbol), reference_point_label,
//                                     '{0:.1f}{1}'.format(float(a[idx]) / quotient, symbol)])

//             else:
//                 xticks_values = [w[idx]]

//                 # only if upstream region is set, add a x tick
//                 if b[idx] > 0:
//                     xticks_values.append(b[idx])
//                     xtickslabel.append('{0:.1f}{1}'.format(-(float(b[idx]) / quotient), symbol))

//                 xtickslabel.append(start_label)

//                 if c[idx] > 0:
//                     xticks_values.append(b[idx] + c[idx])
//                     xtickslabel.append("")

//                 if d[idx] > 0:
//                     xticks_values.append(b[idx] + c[idx] + m[idx])
//                     xtickslabel.append("")

//                 xticks_values.append(b[idx] + c[idx] + m[idx] + d[idx])
//                 xtickslabel.append(end_label)

//                 if a[idx] > 0:
//                     xticks_values.append(b[idx] + c[idx] + m[idx] + d[idx] + a[idx])
//                     xtickslabel.append('{0:.1f}{1}'.format(float(a[idx]) / quotient, symbol))

//                 last = 0
//                 if len(xticks):
//                     last = xticks[-1]
//                 xticks.extend([last + (k / w[idx]) for k in xticks_values])
//         x_axis = np.arange(xticks[-1]) + 1
//         labs = []
//         for x_value in x_axis:
//             if x_value in xticks and xtickslabel[xticks.index(x_value)]:
//                 labs.append(xtickslabel[xticks.index(x_value)])
//             elif x_value in xticks:
//                 labs.append("tick")
//             else:
//                 labs.append("")

//         with open(file_handle, 'w') as fh:
//             # write labels
//             fh.write("bin labels\t\t{}\n".format("\t".join(labs)))
//             fh.write('bins\t\t{}\n'.format("\t".join([str(x) for x in x_axis])))

//             for sample_idx in range(self.matrix.get_num_samples()):
//                 for group_idx in range(self.matrix.get_num_groups()):
//                     sub_matrix = self.matrix.get_matrix(group_idx, sample_idx)
//                     values = [str(x) for x in np.ma.__getattribute__(averagetype)(sub_matrix['matrix'], axis=0)]
//                     fh.write("{}\t{}\t{}\n".format(sub_matrix['sample'], sub_matrix['group'], "\t".join(values)))

//     def save_matrix_values(self, file_name):
//         # print a header telling the group names and their length
//         fh = open(file_name, 'wb')
//         info = []
//         groups_len = np.diff(self.matrix.group_boundaries)
//         for i in range(len(self.matrix.group_labels)):
//             info.append("{}:{}".format(self.matrix.group_labels[i],
//                                        groups_len[i]))
//         fh.write(toBytes("#{}\n".format("\t".join(info))))
//         # add to header the x axis values
//         fh.write(toBytes("#downstream:{}\tupstream:{}\tbody:{}\tbin size:{}\tunscaled 5 prime:{}\tunscaled 3 prime:{}\n".format(
//                  self.parameters['downstream'],
//                  self.parameters['upstream'],
//                  self.parameters['body'],
//                  self.parameters['bin size'],
//                  self.parameters.get('unscaled 5 prime', 0),
//                  self.parameters.get('unscaled 3 prime', 0))))
//         sample_len = np.diff(self.matrix.sample_boundaries)
//         for i in range(len(self.matrix.sample_labels)):
//             info.extend([self.matrix.sample_labels[i]] * sample_len[i])
//         fh.write(toBytes("{}\n".format("\t".join(info))))

//         fh.close()
//         # reopen again using append mode
//         fh = open(file_name, 'ab')
//         np.savetxt(fh, self.matrix.matrix, fmt="%.4g", delimiter="\t")
//         fh.close()

//     def save_BED(self, file_handle):
//         boundaries = np.array(self.matrix.group_boundaries)
//         # Add a header
//         file_handle.write("#chrom\tstart\tend\tname\tscore\tstrand\tthickStart\tthickEnd\titemRGB\tblockCount\tblockSizes\tblockStart\tdeepTools_group")
//         if self.matrix.silhouette is not None:
//             file_handle.write("\tsilhouette")
//         file_handle.write("\n")
//         for idx, region in enumerate(self.matrix.regions):
//             # the label id corresponds to the last boundary
//             # that is smaller than the region index.
//             # for example for a boundary array = [0, 10, 20]
//             # and labels ['a', 'b', 'c'],
//             # for index 5, the label is 'a', for
//             # index 10, the label is 'b' etc
//             label_idx = np.flatnonzero(boundaries <= idx)[-1]
//             starts = ["{0}".format(x[0]) for x in region[1]]
//             ends = ["{0}".format(x[1]) for x in region[1]]
//             starts = ",".join(starts)
//             ends = ",".join(ends)
//             file_handle.write(
//                 '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{1}\t{2}\t0'.format(
//                     region[0],
//                     region[1][0][0],
//                     region[1][-1][1],
//                     region[2],
//                     region[5],
//                     region[4]))
//             file_handle.write(
//                 '\t{0}\t{1}\t{2}\t{3}'.format(
//                     len(region[1]),
//                     ",".join([str(int(y) - int(x)) for x, y in region[1]]),
//                     ",".join([str(int(x) - int(starts[0])) for x, y in region[1]]),
//                     self.matrix.group_labels[label_idx]))
//             if self.matrix.silhouette is not None:
//                 file_handle.write("\t{}".format(self.matrix.silhouette[idx]))
//             file_handle.write("\n")
//         file_handle.close()

//     @staticmethod
//     def matrix_avg(matrix, avgType='mean'):
//         matrix = np.ma.masked_invalid(matrix)
//         return np.ma.__getattribute__(avgType)(matrix, axis=0)

//     def get_individual_matrices(self, matrix):
//         """In case multiple matrices are saved one after the other
//         this method splits them appart.
//         Returns a list containing the matrices
//         """
//         num_cols = matrix.shape[1]
//         num_ind_cols = self.get_num_individual_matrix_cols()
//         matrices_list = []
//         for i in range(0, num_cols, num_ind_cols):
//             if i + num_ind_cols > num_cols:
//                 break
//             matrices_list.append(matrix[:, i:i + num_ind_cols])
//         return matrices_list

//     def get_num_individual_matrix_cols(self):
//         """
//         returns the number of columns  that
//         each matrix should have. This is done because
//         the final matrix that is plotted can be composed
//         of smaller matrices that are merged one after
//         the other.
//         """
//         matrixCols = ((self.parameters['downstream'] + self.parameters['upstream'] + self.parameters['body'] + self.parameters['unscaled 5 prime'] + self.parameters['unscaled 3 prime']) //
//                       self.parameters['bin size'])

//         return matrixCols

// def computeSilhouetteScore(d, idx, labels):
//     """
//     Given a square distance matrix with NaN diagonals, compute the silhouette score
//     of a given row (idx). Each row should have an associated label (labels).
//     """
//     keep = ~np.isnan(d[idx, ])
//     foo = np.bincount(labels[keep], weights=d[idx, ][keep])
//     groupSizes = np.bincount(labels[keep])
//     intraIdx = labels[idx]
//     if groupSizes[intraIdx] == 1:
//         return 0
//     intra = foo[labels[idx]] / groupSizes[intraIdx]
//     interMask = np.arange(len(foo))[np.arange(len(foo)) != labels[idx]]
//     inter = np.min(foo[interMask] / groupSizes[interMask])
//     return (inter - intra) / max(inter, intra)

/// Reference point for heatmap alignment.
enum RefPoint {
    /// Transcription Start Site.
    TSS,
    /// Transcription End Site.
    TES,
    /// Center of the region.
    Center,
}

/// Method for averaging values within a bin.
enum BinAvgType {
    /// Mean value.
    Mean,
    /// Median value.
    Median,
    /// Maximum value.
    Max,
}

/// Arguments for configuring the Heatmapper.
struct HeatmapperArgs {
    /// The bin size to use for the heatmap.
    bin_size: u32,
    /// The number of bins to use upstream of the region.
    upstream: u32,
    /// The number of bins to use downstream of the region.
    downstream: u32,
    /// The number of bins to use in the body of the region.
    body: u32,
    /// The number of bins to use in the unscaled 5 prime region.
    unscaled_5_prime: u32,
    /// The number of bins to use in the unscaled 3 prime region.
    unscaled_3_prime: u32,
    /// The type of average to use when computing the bin values.
    bin_avg_type: BinAvgType,
    /// Whether to treat missing data as zero.
    missing_data_as_zero: bool,
    /// The minimum threshold for a bin value.
    min_threshold: f32,
    /// The maximum threshold for a bin value.
    max_threshold: f32,
    /// The scale to apply to the bin values.
    scale: f32,
    /// The reference point to use when computing the heatmap.
    ref_point: RefPoint,
    /// Whether to treat NaN values after the end of the region as zero.
    nan_after_end: bool,
    /// Whether to print debug information.
    debug: bool,
    /// Whether to print verbose information.
    verbose: bool,
    /// Whether to print quiet information.
    quiet: bool,
}

/// Type of score file (BigWig or BAM).
enum ScoreFileType {
    BigWig,
    Bam,
}

impl FromStr for ScoreFileType {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "bigwig" => Ok(Self::BigWig),
            "bam" => Ok(Self::Bam),
            _ => Err(format!("Invalid score file type: {}", s)),
        }
    }
}

/// Main struct for generating heatmaps from genomic data.
struct Heatmapper {
    score_file: PathBuf,
    score_file_type: ScoreFileType,
    args: HeatmapperArgs,
    regions: Vec<Region>,
    n_regions: usize,
    n_bins: usize,
}

impl Heatmapper {
    /// Creates a new `Heatmapper` instance.
    fn new(score_file: PathBuf, args: HeatmapperArgs, regions: Vec<Region>) -> Self {
        let n_regions = regions.len();
        let n_bins = ((args.upstream + args.downstream + args.body) / args.bin_size) as usize;

        Self {
            score_file: score_file.clone(),
            score_file_type: ScoreFileType::from_str(
                score_file.extension().unwrap().to_str().unwrap(),
            )
            .unwrap(),
            args,
            regions,
            n_regions,
            n_bins,
        }
    }

    /// Reads regions from a BED file.
    fn read_bed(bed: PathBuf) -> Result<Vec<Region>> {
        let mut regions = Vec::new();
        let file = std::fs::File::open(bed)?;
        let r = std::io::BufReader::new(file);
        let mut bed_reader = bed::Reader::new(r);

        for record in bed_reader.records() {
            let record: bed::Record<3> = record?;
            let chrom = record.reference_sequence_name().to_string();
            let start = record.start_position();
            let end = record.end_position();

            regions.push(Region::new(chrom, start..=end));
        }

        Ok(regions)
    }

    /// Creates a `Heatmapper` from a BAM file and a BED file.
    pub fn from_bam(bam: PathBuf, bed: PathBuf, args: HeatmapperArgs) -> Result<Self> {
        // Construct the regions from the BED file
        let regions = Self::read_bed(bed)?;
        let heatmapper = Self::new(bam, args, regions);
        Ok(heatmapper)
    }

    /// Creates a `Heatmapper` from a BigWig file and a BED file.
    pub fn from_bigwig(bigwig: PathBuf, bed: PathBuf, args: HeatmapperArgs) -> Result<Self> {
        // Construct the regions from the BED file
        let regions = Self::read_bed(bed)?;
        let heatmapper = Self::new(bigwig, args, regions);
        Ok(heatmapper)
    }

    /// Computes the submatrix from the BAM file.
    pub fn submatrix_from_bam(&mut self) -> Result<()> {
        let bam = bam::io::indexed_reader::Builder::default().build_from_path(&self.score_file)?;

        // Determine the chunks over which to iterate
        let bam_stats = BamStats::new(self.score_file.clone())?;
        let genomic_chunks = bam_stats.genome_chunks(self.args.bin_size as u64)?;

        // Determine which chunks contain the features of interest
        let genomic_chunk_lapper = regions_to_lapper(genomic_chunks)?;
        let feature_lapper = regions_to_lapper(self.regions.clone())?;

        let genomic_intervals = genomic_chunk_lapper
            .keys()
            .into_iter()
            .map(|chromosome| {
                if feature_lapper.contains_key(chromosome) {
                    let feature_intervals = feature_lapper.get(chromosome).unwrap();
                    let genomic_intervals = genomic_chunk_lapper.get(chromosome).unwrap();

                    // Filter out the genomic intervals that don't contain any features
                    let genomic_intervals_filtered: Vec<Region> = genomic_intervals
                        .intervals
                        .iter()
                        .filter(|iv| feature_intervals.count(iv.start, iv.stop) > 0)
                        .map(|iv| {
                            let start = Position::try_from(iv.start).unwrap();
                            let stop = Position::try_from(iv.stop).unwrap();
                            Region::new(chromosome.clone(), start..=stop)
                        })
                        .collect();
                    Some((chromosome.clone(), genomic_intervals_filtered))
                } else {
                    None
                }
            })
            .filter(|x| x.is_some())
            .map(|x| x.unwrap())
            .collect::<HashMap<String, Vec<Region>>>();

        // Iterate over the genomic intervals and fetch the reads
        genomic_intervals
            .into_par_iter()
            .map(|(chromosome, regions)| {
                let mut reader = bam::io::indexed_reader::Builder::default()
                    .build_from_path(self.score_file.clone())
                    .expect("Error opening BAM file");
                let header = reader.read_header().expect("Error reading BAM header");

                let records = reader
                    .query(&header, &regions)
                    .expect("Error querying BAM file");

                let intervals = records
                    .into_iter()
                    .filter(|r| r.is_ok())
                    .map(|r| r.unwrap())
                    .map(|r| {
                        IntervalMaker::new(
                            r,
                            &header,
                            &chromsizes_refid,
                            &self.filter,
                            self.use_fragment,
                            None,
                        )
                    })
                    .map(|i| i.coords())
                    .filter(|c| c.is_some())
                    .map(|c| c.unwrap())
                    .map(|i| Iv {
                        start: i.0,
                        stop: i.1,
                        val: 1,
                    })
                    .collect::<Vec<Iv>>();
            });

        Ok(())
    }

    // let pileup = genomic_chunks
    //     .into_par_iter()
    //     .progress_with(progress_bar(
    //         n_total_chunks as u64,
    //         "Performing pileup".to_string(),
    //     ))
    //     .map(|region| {
    //         // Open the BAM file
    //         let mut reader = bam::io::indexed_reader::Builder::default()
    //             .build_from_path(self.file_path.clone())
    //             .expect("Error opening BAM file");
    //         // Extract the header
    //         let header = reader.read_header().expect("Error reading BAM header");

    //         // Fetch the reads in the region
    //         let records = reader
    //             .query(&header, &region)
    //             .expect("Error querying BAM file");

    //         // Make intervals from the reads in the region
    //         let intervals = records
    //             .into_iter()
    //             .filter(|r| r.is_ok())
    //             .map(|r| r.unwrap())
    //             .map(|r| {
    //                 IntervalMaker::new(
    //                     r,
    //                     &header,
    //                     &chromsizes_refid,
    //                     &self.filter,
    //                     self.use_fragment,
    //                     None,
    //                 )
    //             })
    //             .map(|i| i.coords())
    //             .filter(|c| c.is_some())
    //             .map(|c| c.unwrap())
    //             .map(|i| Iv {
    //                 start: i.0,
    //                 stop: i.1,
    //                 val: 1,
    //             })
    //             .collect::<Vec<Iv>>();

    //         // Create a lapper from the reads
    //         let mut read_lapper = Lapper::new(intervals);

    //         // Iterate over bins within the given region and count the reads using the lapper
    //         let mut bin_counts: Vec<Iv> = Vec::new();

    //         // Generate the bins for the region (bin_size)
    //         let region_interval = region.interval();
    //         let region_start = region_interval
    //             .start()
    //             .expect("Error getting interval start")
    //             .get();
    //         let region_end = region_interval
    //             .end()
    //             .expect("Error getting interval end")
    //             .get();
}
