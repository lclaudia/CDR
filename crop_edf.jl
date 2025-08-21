using EDF
using EDF: TimestampedAnnotationList, PatientID, RecordingID, SignalHeader,
           Signal, AnnotationsSignal
using Dates
using FilePathsBase

# modified from https://github.com/beacon-biosignals/EDF.jl/blob/main/src/write.jl
# EDF+C is for edf+ files  -  removed
function EDF.write_header(io::IO, file::EDF.File)
    length(file.signals) <= 9999 ||
        error("EDF does not allow files with more than 9999 signals")
    expected_bytes_written = EDF.BYTES_PER_FILE_HEADER +
                             EDF.BYTES_PER_SIGNAL_HEADER * length(file.signals)
    bytes_written = 0
    bytes_written += EDF.edf_write(io, file.header.version, 8)
    bytes_written += EDF.edf_write(io, file.header.patient, 80)
    bytes_written += EDF.edf_write(io, file.header.recording, 80)
    bytes_written += EDF.edf_write(io, file.header.start, 16)
    bytes_written += EDF.edf_write(io, expected_bytes_written, 8)
    ####
    bytes_written += EDF.edf_write(io, file.header.is_contiguous ? "     " : "EDF+D", 44)
    # bytes_written += EDF.edf_write(io, file.header.is_contiguous ? "EDF+C" : "EDF+D", 44)
    ####
    bytes_written += EDF.edf_write(io, file.header.record_count, 8)
    bytes_written += EDF.edf_write(io, file.header.seconds_per_record, 8)
    bytes_written += EDF.edf_write(io, length(file.signals), 4)
    signal_headers = EDF.SignalHeader.(file.signals)
    for (field_name, byte_limit) in EDF.SIGNAL_HEADER_FIELDS
        for signal_header in signal_headers
            field = getfield(signal_header, field_name)
            bytes_written += EDF.edf_write(io, field, byte_limit)
        end
    end
    bytes_written += EDF.edf_write(io, ' ', 32 * length(file.signals))
    @assert bytes_written == expected_bytes_written
    return bytes_written
end



EDF_file = "long_EDF_file.edf";               # use her the file name of the original recording
                                              # this patient has 78 channels
CH=2:79;                                      # set here the channel numbers from the original edf file  
SR = 500;                                     # set the sampling rate
AnfEnde = [18589000 18709000];                # start end end of the data segment to be extracted

edf=EDF.read(EDF_file);                       # read the file 

daten=edf.signals;
hdr=edf.header;

edf_A_ori = length(daten[CH[1]].samples);
edf_B_ori = daten[CH[1]].header.samples_per_record;
edf_C_ori = hdr.record_count;

if edf_B_ori*edf_C_ori != edf_A_ori           # sometimes this needs to be corrected
   edf_B_ori=Int(edf_A_ori/edf_C_ori)
end

AE=AnfEnde[1]:AnfEnde[2]; L_AE=length(AE);

edf_A = L_AE;
edf_B = edf_B_ori;
edf_C = Int(ceil(edf_A/edf_B));
edf_A = edf_B*edf_C;
edf_D = edf_B/SR;

SNR=15;

X=fill(0,length(CH)*2,edf_A);
for ch=1:length(CH)
                                              # data segment
    X[ch,1:L_AE] = daten[CH[ch]].samples[AE];
                                              # remove mean
    X[ch,1:L_AE] .-= Int(round(mean(X[ch,1:L_AE])));
                                              # add noise and convert to integer
    X[length(CH)+ch,1:L_AE] = Int.(round.(add_noise(X[ch,1:L_AE],SNR)));
end

header = EDF.FileHeader(
       hdr.version,
       #hdr.patient,
       "patient 1",                           # new patient name
       #hdr.recording,                        # new data
       "Date:01.01.10 Time:00.00.00",
       #hdr.start,
       DateTime("2010-01-01T00:00:00"),
       hdr.is_contiguous,
       #hdr.record_count,
       edf_C,                                 # number of segments for each channel
       #hdr.seconds_per_record
       edf_D                                  # seconds_per_record
       );

signals = Array{Union{EDF.AnnotationsSignal, EDF.Signal{Int16}},1}(undef,length(CH)*2);
for i in 1:length(CH)
    signal_header = EDF.SignalHeader(
        daten[CH[i]].header.label,
        daten[CH[i]].header.transducer_type,                     
        daten[CH[i]].header.physical_dimension, 
        daten[CH[i]].header.physical_minimum,
        daten[CH[i]].header.physical_maximum,
        daten[CH[i]].header.digital_minimum,
        daten[CH[i]].header.digital_maximum,
        daten[CH[i]].header.prefilter,
        edf_B
        );
    signals[i] = EDF.Signal{Int16}(signal_header,X[i,:]);
                                              # the 78 channels 

    signal_header = EDF.SignalHeader(
        @sprintf("%s_%02ddB",signal_header.label,SNR),
        daten[CH[i]].header.transducer_type,                     
        daten[CH[i]].header.physical_dimension, 
        daten[CH[i]].header.physical_minimum,
        daten[CH[i]].header.physical_maximum,
        daten[CH[i]].header.digital_minimum,
        daten[CH[i]].header.digital_maximum,
        daten[CH[i]].header.prefilter,
        edf_B
        );
    signals[length(CH)+i] = EDF.Signal{Int16}(signal_header,X[length(CH)+i,:]);
                                              # the 78 channels with added noise
end

EDF_file_out = "patient1_S05__01_03.edf";

open(EDF_file_out, "w") do io
    edf_file = EDF.File(io, header, signals);
    EDF.write(io, edf_file);
end
#the do-block form of open can be used to automatically close the file even in the case of exceptions
