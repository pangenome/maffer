"""
mafValidator

Script to validate Multiple Alignment Format (maf) files.

Python 3 porting of mafValidator (https://github.com/dentearl/mafTools/blob/master/mafValidator/src/mafValidator.py)
from mafTools (https://github.com/dentearl/mafTools)
"""

import os
import numpy as np
import re
import sys
import argparse

g_version = '0.1.1 Sept 2020'


class ValidatorError(Exception):
    pass


class SourceLengthError(ValidatorError):
    pass


class SpeciesFieldError(ValidatorError):
    pass


class MissingAlignmentBlockLineError(ValidatorError):
    pass


class KeyValuePairError(ValidatorError):
    pass


class AlignmentLengthError(ValidatorError):
    pass


class FieldNumberError(ValidatorError):
    pass


class FooterError(ValidatorError):
    pass


class StrandCharacterError(ValidatorError):
    pass


class StartFieldError(ValidatorError):
    pass


class SourceSizeFieldError(ValidatorError):
    pass


class OutOfRangeError(ValidatorError):
    pass


class HeaderError(ValidatorError):
    pass


class ILineFormatError(ValidatorError):
    pass


class ELineFormatError(ValidatorError):
    pass


class QLineFormatError(ValidatorError):
    pass


class DuplicateColumnError(ValidatorError):
    pass


class SequenceConsistencyError(ValidatorError):
    pass


class EmptyInputError(ValidatorError):
    pass


class GenericValidationOptions:
    # used by other python modules when calling mafValidator from unit tests
    def __init__(self):
        self.lookForDuplicateColumns = False
        self.testChromNames = False
        self.validateSequence = True


def init_options(parser):
    parser.add_argument('--maf', dest='filename',
                        help='path to maf file to validate.')
    parser.add_argument('--testChromNames', dest='testChromNames',
                        action='store_true',
                        default=False,
                        help='Test that all species fields contain chrom name, i.e.: s hg19.chr1')
    parser.add_argument('--ignoreDuplicateColumns', dest='lookForDuplicateColumns',
                        default=True, action='store_false',
                        help=('Turn off the checks for duplicate columns, may be useful for pairwise-only '
                              'alignments. default=duplicate checking is on.'))
    parser.add_argument('--validateSequence', dest='validateSequence',
                        default=False, action='store_true',
                        help=('Turn on checks to make sure all sequence fields are '
                              'consistent. Slows things down considerably. Note that selecting this option '
                              'implicitly sets --ignoreDuplicateColumns'))
    parser.add_argument('--version', dest='isVersion', action='store_true', default=False,
                        help='Print version number and exit.')


def check_options(args, parser):
    if args.isVersion:
        print('mafValidator.py, version %s' % g_version)
        sys.exit(0)
    if args.filename is None:
        parser.error('specify --maf')
    if not os.path.exists(args.filename):
        parser.error('--maf %s does not exist.' % args.filename)
    if args.validateSequence:
        args.lookForDuplicateColumns = False


def validateILine(lineno, line, prevline, filename):
    """ Checks all lines that start with 'i' and raises an exepction if
    the line is malformed.
    i lines are made up of five fields after the 'i':
    From http://genome.ucsc.edu/FAQ/FAQformat#format5 :
    src -- The name of the source sequence for the alignment.
           Should be the same as the 's' line immediately above this line.
    leftStatus -- A character that specifies the relationship between the
                  sequence in this block and the sequence that appears in
                  the previous block.
    leftCount -- Usually the number of bases in the aligning species between
                 the start of this alignment and the end of the previous one.
    rightStatus -- A character that specifies the relationship between the
                   sequence in this block and the sequence that appears in
                   the subsequent block.
    rightCount -- Usually the number of bases in the aligning species between
                  the end of this alignment and the start of the next one.
    The status characters can be one of the following values:
    C -- the sequence before or after is contiguous with this block.
    I -- there are bases between the bases in this block and the one before or after it.
    N -- this is the first sequence from this src chrom or scaffold.
    n -- this is the first sequence from this src chrom or scaffold but it is
         bridged by another alignment from a different chrom or scaffold.
    M -- there is missing data before or after this block (Ns in the sequence).
    T -- the sequence in this block has been used before in a previous block (likely a tandem duplication)
    """
    d = line.split()
    p = prevline.split()
    if len(d) != 6:
        raise ILineFormatError('maf %s contains an "i" line that has too many fields on line number %d: '
                               '%s' % (filename, lineno, line))
    if p[0] != 's':
        raise ILineFormatError('maf %s contains an "i" line that does not follow an "s" line on line number %d: '
                               '%s' % (filename, lineno, line))
    for i in [3, 5]:
        try:
            n = int(d[i])
        except ValueError:
            raise ILineFormatError('maf %s contains an "i" line that has non integer Count "%s" on line number %d: '
                                   '%s' % (filename, d[i], lineno, line))
        if int(d[i]) < 0:
            raise ILineFormatError('maf %s contains an "i" line that has negative Count "%s" on line number %d: '
                                   '%s' % (filename, d[i], lineno, line))
    for i in [2, 4]:
        if d[i] not in ['C', 'I', 'N', 'n', 'M', 'T']:
            raise ILineFormatError('maf %s contains an "i" line with an invalid Status "%s" on line number %d: '
                                   '%s' % (filename, d[i], lineno, line))
        if d[i] == 'I' and int(d[i + 1]) < 1:
            raise ILineFormatError('maf %s contains an "i" line with an invalid Count "%s" on line number %d: '
                                   '%s' % (filename, d[i + 1], lineno, line))
    if p[1] != d[1]:
        raise ILineFormatError('maf %s contains an "i" line with a different src value "%s" than on the previous '
                               '"s" line "%s" on line number %d: '
                               '%s' % (filename, d[1], p[1], lineno, line))


def validateELine(lineno, line, filename):
    """ Checks all lines that start with 'e' and raises an exepction if
    the line is malformed.
    e lines are made up of six fields after the 'e':
    From http://genome.ucsc.edu/FAQ/FAQformat#format5 :

    src -- The name of one of the source sequences for the alignment.
    start -- The start of the non-aligning region in the source sequence.
             This is a zero-based number. If the strand field is '-' then
             this is the start relative to the reverse-complemented source
             sequence (see Coordinate Transforms).
    size -- The size in base pairs of the non-aligning region in the source sequence.
    strand -- Either '+' or '-'. If '-', then the alignment is to the reverse-complemented source.
    srcSize -- The size of the entire source sequence, not just the parts
               involved in the alignment. alignment and any insertions (dashes) as well.
    status -- A character that specifies the relationship between the non-aligning
              sequence in this block and the sequence that appears in the previous and subsequent blocks.
    The status character can be one of the following values:
    C -- the sequence before and after is contiguous implying that this region
         was either deleted in the source or inserted in the reference sequence.
         The browser draws a single line or a '-' in base mode in these blocks.
    I -- there are non-aligning bases in the source species between chained
         alignment blocks before and after this block. The browser shows a double line or '=' in base mode.
    M -- there are non-aligning bases in the source and more than 90% of them
         are Ns in the source. The browser shows a pale yellow bar.
    n -- there are non-aligning bases in the source and the next aligning block
         starts in a new chromosome or scaffold that is bridged by a chain between
         still other blocks. The browser shows either a single line or a double
         line based on how many bases are in the gap between the bridging alignments.
    """
    d = line.split()
    if len(d) != 7:
        raise ELineFormatError('maf %s contains an "e" line that has too many fields on line number %d: '
                               '%s' % (filename, lineno, line))
    for i in [2, 3, 5]:
        try:
            n = int(d[i])
        except ValueError:
            raise ELineFormatError('maf %s contains an "e" line that has non integer Count "%s" on line number %d: '
                                   '%s' % (filename, d[i], lineno, line))
        if int(d[i]) < 0:
            raise ELineFormatError('maf %s contains an "e" line that has negative Count "%s" on line number %d: '
                                   '%s' % (filename, d[i], lineno, line))
    if d[4] not in ['+', '-']:
        raise ELineFormatError('maf %s contains an "e" line with an invalid strand "%s" on line number %d: '
                               '%s' % (filename, d[4], lineno, line))
    if d[6] not in ['C', 'I', 'M', 'n']:
        raise ELineFormatError('maf %s contains an "e" line with an invalid Status "%s" on line number %d: '
                               '%s' % (filename, d[2], lineno, line))


def validateQLine(lineno, line, prevline, filename):
    """ Checks all lines that start with 'q' and raises an exepction if
    the line is malformed.
    q lines are made up of two fields after the 'q':
    From http://genome.ucsc.edu/FAQ/FAQformat#format5 :
    The 'q' lines contain a compressed version of the actual raw quality data,
    representing the quality of each aligned base for the species with a single
    character of 0-9 or F. The following fields are defined by position rather
    than name=value pairs:
    src -- The name of the source sequence for the alignment. Should be the same
           as the 's' line immediately preceding this line.
    value -- A MAF quality value corresponding to the aligning nucleotide acid
             in the preceding 's' line. Insertions (dashes) in the preceding 's'
             line are represented by dashes in the 'q' line as well. The quality
             value can be 'F' (finished sequence) or a number derived from the
             actual quality scores (which range from 0-97) or the manually
             assigned score of 98. These numeric values are calculated as:
             MAF quality value = min( floor(actual quality value/5), 9 )
    This results in the following mapping:
    MAF quality value | Raw quality score range | Quality level
    0-8 |  0-44 | Low
      9 | 45-97 | High
      0 |    98 | Manually assigned
      F |    99 | Finished

    """
    d = line.split()
    p = prevline.split()
    if len(d) != 3:
        raise QLineFormatError('maf %s contains an "q" line that has too many fields on line number %d: '
                               '%s' % (filename, lineno, line))
    if p[0] != 's':
        raise QLineFormatError('maf %s contains an "q" line that does not follow an "s" line on line number %d: '
                               '%s' % (filename, lineno, line))
    for c in d[2]:
        if c not in map(str, range(0, 10)) + ['F', '-']:
            raise QLineFormatError('maf %s contains an "q" line with an invalid character "%s" on line number %d: '
                                   '%s' % (filename, c, lineno, line))
    if p[1] != d[1]:
        raise QLineFormatError('maf %s contains an "q" line with a different src value "%s" than on the previous '
                               '"s" line "%s" on line number %d: '
                               '%s' % (filename, d[1], p[1], lineno, line))


def validateAlignmentLine(lineno, line, filename):
    """ Checks all lines that start with 'a' and raises an exception if
    the line is malformed.
    """
    validateKeyValuePairLine(lineno, line, filename)


def validateKeyValuePairLine(lineno, line, filename):
    d = line.split()
    for i in range(1, len(d)):
        if len(d[i].split('=')) != 2:
            raise KeyValuePairError('maf %s has a line that does not contain '
                                    'good key-value pairs on line number %d: %s'
                                    % (filename, lineno, line))


def validateSeqLine(namePat, options, lineno, line, filename, sequenceColumnDict, sequenceFieldDict):
    data = line.split()
    if len(data) != 7:
        raise FieldNumberError('maf %s has incorrect number of fields on line number %d: %s'
                               % (filename, lineno, line))
    if data[4] not in ('-', '+'):
        raise StrandCharacterError('maf %s has unexpected character in strand field "%s" on line number %d: %s'
                                   % (filename, data[4], lineno, line))
    if data[4] == '+':
        strand = 1
    else:
        strand = -1
    if int(data[3]) != len(data[6].replace('-', '')):
        raise AlignmentLengthError(
            'maf %s sequence length field (%d) contradicts alignment field (non-gapped length %d) on line number %d: %s'
            % (filename, int(data[3]), len(data[6].replace('-', '')), lineno, line))
    if int(data[2]) < 0:
        raise StartFieldError('maf %s has bad start field on line number %d: %s'
                              % (filename, lineno, line))
    if int(data[5]) < 0:
        raise SourceSizeFieldError('maf %s has bad source size field on line number %d: %s'
                                   % (filename, lineno, line))
    if int(data[2]) + int(data[3]) > int(data[5]):
        raise OutOfRangeError('maf %s out of range sequence on line number %d: %s'
                              % (filename, lineno, line))
    if options.lookForDuplicateColumns:
        checkForDuplicateColumns(data, sequenceColumnDict, filename, lineno, line)
    if options.validateSequence:
        checkForSequenceInconsistencies(data, sequenceFieldDict, filename, lineno, line)
    if options.testChromNames:
        m = re.match(namePat, data[1])
        if m is None:
            raise SpeciesFieldError('maf %s has name (source) field without ".chr" suffix: "%s" on line number %d: %s'
                                    % (filename, data[1], lineno, line))
        return m.group(1), m.group(2), data[5], len(data[6])
    return data[1], None, data[5], len(data[6])


def reverseComplement(s):
    s = s[::-1]
    s = s.replace('A', '1')
    s = s.replace('a', '1')
    s = s.replace('T', 'A')
    s = s.replace('t', 'A')
    s = s.replace('1', 'T')
    s = s.replace('G', '2')
    s = s.replace('g', '2')
    s = s.replace('C', 'G')
    s = s.replace('c', 'G')
    s = s.replace('2', 'C')
    return s


def checkForSequenceInconsistencies(data, sfd, filename, lineno, line):
    """ sfd = sequenceFieldDict, a dictionary keyed on sequence field names and valued with
    numpy arrays (dtype chararray) that stores all representations of the sequence in the maf.
    If the sequence should change over the course of the file, throws an error.
    """
    name = data[1]
    start, length = map(int, data[2:4])
    strand = data[4]
    total_src_length = int(data[5])
    if name not in sfd:
        sfd[name] = np.zeros(int(total_src_length), dtype=np.str)
        sfd[name][:] = 'N'

    if strand == '+':
        stop = start + length
        seq_str = data[6]
    else:
        start, stop = total_src_length - (start + length), total_src_length - start
        seq_str = reverseComplement(data[6])
    print(data)
    print(name, start, stop, seq_str)
    # remove gaps, switch to uppercase
    seq_str = seq_str.upper().replace('-', '')
    seq = np.zeros(length, dtype=np.str)
    seq[:] = list(seq_str)
    stored = np.ma.array(sfd[name][start:stop], mask=sfd[name][start:stop] == 'N')
    if not (stored == seq).all():
        if not (stored == seq).mask.all():
            if stored is not np.ma.masked:
                raise SequenceConsistencyError('maf %s has inconsistent sequence, discovered on line number %d. '
                                               '\nFirst: %s\nNow  : %s\n%s'
                                               % (filename, lineno, stored, seq_str, line))
    sfd[name][start:stop] = seq
    print(sfd[name])
    print('-------------')

def checkForDuplicateColumns(data, scd, filename, lineno, line):
    """ scd = sequenceColumnDict, a dictionary keyed on sequence field names and valued with
    numpy arrays (dtype boolean) that indicates whether or not a column has already appeared in
    the alignment.
    """
    # data[1] sequence name
    # data[2] start position
    # data[3] seq length
    # data[4] strand
    # data[5] total source length
    name = data[1]
    start, length = map(int, data[2:4])
    strand = data[4]
    total_src_length = int(data[5])
    if name not in scd:
        scd[name] = np.zeros(int(total_src_length), dtype=np.bool)
    if strand == '+':
        stop = start + length
    else:
        start, stop = total_src_length - (start + length), total_src_length - start
    if scd[name][start:stop].any():
        raise DuplicateColumnError('maf %s has duplicate columns, second instance discovered on line number %d: %s'
                                   % (filename, lineno, line))
    scd[name][start:stop] = True


def validate_header(f, filename):
    """ tests the first line of the maf file to make sure it is valid.
    valid starts are either "track ..." or "##maf..."
    """
    header = f.readline()
    lineno = 1
    if header.startswith('track'):
        try:
            validateKeyValuePairLine(1, header, filename)
        except KeyValuePairError:
            raise HeaderError('maf %s has bad header, fails key value '
                              'pair tests, %s on linenumber %d' % (filename, header, lineno))
        header = f.readline()
        lineno += 1
    if not header.startswith('##maf'):
        raise HeaderError('maf %s has bad header, fails to start with `##\': %s'
                          % (filename, header))
    try:
        validateKeyValuePairLine(lineno, header, filename)
    except KeyValuePairError:
        raise HeaderError('maf %s has bad header, fails key value '
                          'pair tests, %s on linenumber %d' % (filename, header, lineno))
    data = header.split()
    version = False
    for d in data[1:]:
        if d.startswith('=') or d == '=' or d.endswith('='):
            raise HeaderError('maf %s has bad header, there should be '
                              'no whitespace surrounding "=": %s'
                              % (filename, header))
        try:
            k, v = d.split('=')
        except ValueError:
            raise HeaderError('maf %s has bad header, there should be '
                              'no whitespace surrounding "=": %s'
                              % (filename, header))
        if k == 'version':
            version = True
    if not version:
        raise HeaderError('maf %s has bad header, no version information: %s'
                          % (filename, header))
    return lineno


def validate_maf(filename, options):
    """ returns true on valid maf file
    a completely empty file should be considered invalid.
    """
    if os.path.getsize(filename) == 0:
        # empty files are considered invalid
        raise EmptyInputError('maf %s is completely empty.' % filename)
    name_regex = r'(.+?)\.(chr.+)'
    name_pat = re.compile(name_regex)
    sequence_column_dict = {}
    sequence_field_dict = {}

    with open(filename, 'rt') as f:
        header_processed_lines = validate_header(f, filename)
        sources = {}
        prev_line_was_alignment_block = False
        alignment_field_length = None
        prev_line = ''
        for lineno, line in enumerate(f, header_processed_lines + 1):
            line = line.strip()
            if line.startswith('#'):
                alignment_field_length = None
                prev_line = line
                continue
            elif line.startswith('a'):
                validateAlignmentLine(lineno, line, filename)
                prev_line_was_alignment_block = True
                alignment_field_length = None
            elif line.startswith('s'):
                if not prev_line_was_alignment_block:
                    raise MissingAlignmentBlockLineError('maf %s has a sequence line that was not preceded '
                                                         'by an alignment line on line number %d: %s'
                                                         % (filename, lineno, line))
                species, chrom, length, al_field_len = validateSeqLine(name_pat, options, lineno, line,
                                                                       filename, sequence_column_dict,
                                                                       sequence_field_dict)
                if alignment_field_length is None:
                    alignment_field_length = al_field_len
                else:
                    if alignment_field_length != al_field_len:
                        raise AlignmentLengthError('maf %s has a sequence line with an alignment field of different '
                                                   'length than the other sequences in the block on line number %d: %s'
                                                   % (filename, lineno, line))
                if (species, chrom) in sources:
                    if sources[(species, chrom)] != length:
                        if chrom is None:
                            raise SourceLengthError('maf %s has different source lengths for '
                                                    'species %s (read %s but was previously %s) on line number %d: %s'
                                                    % (
                                                        filename, species, sources[(species, chrom)], length, lineno,
                                                        line))
                        else:
                            raise SourceLengthError('maf %s has different source lengths for '
                                                    'species %s chrom %s (read %s but was previously %s) on line '
                                                    'number %d: %s '
                                                    % (
                                                        filename, species, chrom, sources[(species, chrom)], length,
                                                        lineno,
                                                        line
                                                    ))
                else:
                    sources[(species, chrom)] = length
            elif line.startswith('i'):
                if not prev_line_was_alignment_block:
                    raise MissingAlignmentBlockLineError('maf %s has an "i" line that was not preceded '
                                                         'by an alignment line on line number %d: %s'
                                                         % (filename, lineno, line))
                validateILine(lineno, line, prev_line, filename)
            elif line.startswith('e'):
                if not prev_line_was_alignment_block:
                    raise MissingAlignmentBlockLineError('maf %s has a "e" line that was not preceded '
                                                         'by an alignment line on line number %d: %s'
                                                         % (filename, lineno, line))
                validateELine(lineno, line, filename)
            elif line.startswith('q'):
                if not prev_line_was_alignment_block:
                    raise MissingAlignmentBlockLineError('maf %s has a "q" line that was not preceded '
                                                         'by an alignment line on line number %d: %s'
                                                         % (filename, lineno, line))
                validateQLine(lineno, line, prev_line, filename)
            elif line == '':
                prev_line_was_alignment_block = False
                alignment_field_length = None
            prev_line = line

        if line != '':
            if prev_line_was_alignment_block:
                raise FooterError('maf %s has a bad footer, last alignment block was not closed with a new line.'
                                  % filename)
            else:
                if not line.startswith('#'):
                    raise FooterError('maf %s has a bad footer, should end with a blank line: %s'
                                      % (filename, line))
    return True


def main():
    parser = argparse.ArgumentParser()

    init_options(parser)
    args = parser.parse_args()
    check_options(args, parser)

    validate_maf(args.filename, args)


if __name__ == '__main__':
    main()
