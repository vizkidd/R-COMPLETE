#include "omp.h"
#include <arrow/api.h>
#include <arrow/filesystem/localfs.h>
#include <arrow/ipc/writer.h>
// #include <arrow/ipc/type_fwd.h>
#include <arrow/ipc/options.h>
#include <arrow/io/api.h>
#include <arrow/filesystem/filesystem.h>
#include <parquet/arrow/reader.h>
#include <parquet/arrow/writer.h>
#include <boost/lexical_cast.hpp>

class ArrowWrapper
{
private:
    static std::shared_ptr<ArrowWrapper> self;
    std::shared_ptr<arrow::RecordBatchVector> recBth_vec;
    std::shared_ptr<arrow::Field> fasta_field0, fasta_field1, fasta_field2;
    std::shared_ptr<arrow::Schema> fasta_schema;
    std::shared_ptr<arrow::Buffer> buffer;
    std::shared_ptr<arrow::KeyValueMetadata> fasta_metadata;
    // arrow::ipc::RecordBatchStreamReader ipc_stream_reader;
    // arrow::ipc::RecordBatchFileReader ipc_file_reader;
    arrow::ipc::IpcWriteOptions ipc_options;
    std::shared_ptr<std::ostringstream> outputStream;
    ArrowWrapper(/* args */);

public:
    ~ArrowWrapper();
    static std::shared_ptr<ArrowWrapper> GetSelf();
    std::shared_ptr<arrow::Schema> GetFASTASchema(void);
    std::shared_ptr<arrow::KeyValueMetadata> GetFASTAMetadata(void);
    std::vector<std::shared_ptr<arrow::Field>> GetFASTAFields(void);
    arrow::ipc::IpcWriteOptions ArrowWrapper::GetArrowIPCOptions(void);
    static void SetSelf(ArrowWrapper *_obj);
    std::shared_ptr<arrow::ipc::RecordBatchWriter> GetRecordBatchWriter(std::shared_ptr<arrow::io::OutputStream> outStream);
    std::shared_ptr<arrow::RecordBatch> RecordBatchFromString(const std::string_view &str);
    std::shared_ptr<arrow::RecordBatchVector> RecordBatchVectorFromFile(const std::string_view &filename, const int &num_threads);
    int GetColumnCount(const std::string_view &filename, char delim = '\t');
    std::shared_ptr<arrow::RecordBatch> ArrowWrapper::RecordBatchFrom(const FastaSequenceData &input);
    template <typename T>
    std::shared_ptr<arrow::RecordBatch> ArrowWrapper::RecordBatchFrom(const T &str);
    std::shared_ptr<std::vector<std::shared_ptr<std::ostringstream>>> SplitFileIntoStreams(const std::string_view &filename, const int &num_threads, const int &number_of_streams = 1, const char *delim = "");
};

ArrowWrapper::ArrowWrapper()
{

    ipc_options.use_threads = true;
    ipc_options.allow_64bit = true;
    ipc_options.metadata_version = arrow::ipc::MetadataVersion::V5;
    ipc_options.codec = arrow::util::Codec::Create(arrow::Compression::LZ4_FRAME).ValueOrDie();
    ipc_options.unify_dictionaries = true;
    ipc_options.emit_dictionary_deltas = true;
    fasta_metadata = std::make_shared<arrow::KeyValueMetadata>();
    fasta_metadata->Append("format", "FASTA");
    fasta_metadata->Append("type", "nucleotide");

    fasta_field0 = arrow::field("index", arrow::int32());
    fasta_field1 = arrow::field("header", arrow::utf8());
    fasta_field2 = arrow::field("sequence", arrow::utf8());
    fasta_schema = arrow::schema({fasta_field0, fasta_field1, fasta_field2});
    buffer = arrow::Buffer::Wrap("", 0);
    SetSelf(this);
}

ArrowWrapper::~ArrowWrapper()
{
}

static std::shared_ptr<ArrowWrapper> ArrowWrapper::GetSelf()
{
    if (ArrowWrapper::self.get() == nullptr)
    {
        ArrowWrapper *_obj = new ArrowWrapper();
        ArrowWrapper::SetSelf(_obj);
        // BLASTWrapper::self = static_cast<std::shared_ptr<BLASTWrapper>>(_obj);
    }
    return (ArrowWrapper::self);
}

static void ArrowWrapper::SetSelf(ArrowWrapper *_obj)
{
    self = std::shared_ptr<ArrowWrapper>(_obj);
}

std::shared_ptr<arrow::ipc::RecordBatchWriter> ArrowWrapper::GetRecordBatchWriter(std::shared_ptr<arrow::io::OutputStream> outStream)
{
    return arrow::ipc::MakeStreamWriter(outStream.get(), fasta_schema, ipc_options).ValueOrDie();
}

std::shared_ptr<std::vector<std::shared_ptr<std::ostringstream>>> ArrowWrapper::SplitFileIntoStreams(const std::string_view &filename, const int &num_threads, const int &number_of_streams, const char *delim)
{
    // Debug
    std::cout << std::stoi(delim) << std::endl;

    std::shared_ptr<std::vector<std::shared_ptr<std::ostringstream>>> outputStreams = std::make_shared<std::vector<std::shared_ptr<std::ostringstream>>>();
    omp_lock_t outputStreamsLock;

    omp_init_lock(&outputStreamsLock);
#pragma omp parallel num_threads(num_threads) shared(outputStreams)
    {
#pragma omp for nowait schedule(auto)
        for (int i = 0; i < number_of_streams; i++)
        {
            omp_set_lock(&outputStreamsLock);
            outputStreams->push_back(std::make_shared<std::ostringstream>());
            omp_unset_lock(&outputStreamsLock);
        }
#pragma omp barrier
    }

    std::shared_ptr<FILE> file(fopen(filename.data(), "rb"), [](FILE *f)
                               { fclose(f); });
    if (!file)
    {
        std::cerr << "Error: Failed to open file: " << filename << std::endl;
        return nullptr;
    }

    // Get the file size
    fseek(file.get(), 0, SEEK_END);
    long fileSize = ftell(file.get());
    rewind(file.get());

    // Memory map the file
    char *fileDataPtr = static_cast<char *>(mmap(nullptr, fileSize, PROT_READ, MAP_PRIVATE, fileno(file.get()), 0));

    if (fileDataPtr == MAP_FAILED)
    {
        std::cerr << "Error: Failed to map file: " << filename << std::endl;
        fclose(file.get());
        return nullptr;
    }

    std::shared_ptr<char> fileData(fileDataPtr, [fileSize](char *ptr)
                                   { munmap(ptr, fileSize); });

    std::shared_ptr<char> end_of_file = std::shared_ptr<char>(fileData, fileData.get() + fileSize); // fileData.get() + fileSize;
    // Process the FASTA file
    // std::shared_ptr<char> p = fileData;
    // std::shared_ptr<char> seeker(strchr(fileData.get(), '>'));
    char *_seek = strchr(fileData.get(), std::stoi(delim));
    if (_seek == nullptr)
    {
        std::cerr << "Error: Invalid file: Could not seek" << std::to_string(*delim) << ":" << filename << std::endl;
        return nullptr;
    }

    char *p = _seek;
    omp_lock_t pLock;

    omp_init_lock(&pLock);
#pragma omp parallel num_threads(number_of_streams) shared(outputStreams) shared(end_of_file) shared(p)
    {
#pragma omp for schedule(auto) nowait
        for (p = _seek; p < end_of_file.get(); p++)
        {
            omp_set_lock(&pLock);
            if (p != nullptr)
            {
                // Read the header line
                char *entryStart = p;
                char *entryEnd = strchr(p, std::stoi(delim));
                // if (entryStart != nullptr)
                //{
                if (entryEnd == entryStart)
                {
                    // Read once more
                    entryEnd = strchr(p + 1, std::stoi(delim));
                }

                if (entryEnd == nullptr)
                {
                    entryEnd = end_of_file.get();
                }
                if (entryStart == p)
                {
                    std::string_view full_entry(entryStart, entryEnd - entryStart - 1);
                    if (!full_entry.empty())
                    {
                        int omp_thread_id = omp_get_thread_num();
                        omp_set_lock(&outputStreamsLock);
                        std::shared_ptr<std::ostringstream> outputStream = outputStreams->at(omp_thread_id);
                        outputStream->write(full_entry.data(), full_entry.size());
                        omp_unset_lock(&outputStreamsLock);
                    }
                }
                p = entryEnd - 1;
                omp_unset_lock(&pLock);
            }
        }
    }

#pragma omp barrier
    omp_destroy_lock(&pLock);
    omp_destroy_lock(&outputStreamsLock);
    return outputStreams;
}

std::shared_ptr<arrow::RecordBatchVector> ArrowWrapper::RecordBatchVectorFromFile(const std::string_view &filename, const int &num_threads)
{
}

template <>
std::shared_ptr<arrow::RecordBatch> ArrowWrapper::RecordBatchFrom(const FastaSequenceData &input)
{
    std::shared_ptr<arrow::StringBuilder> fasta_h_builder, fasta_seq_builder;
    std::shared_ptr<arrow::Int32Builder> rec_no_builder;
    std::cout << "Processing FastaSequenceData:" << std::endl;
    std::cout << "Header: " << input.header << std::endl;
    std::cout << "Sequence: " << input.seq << std::endl;
    fasta_seq_builder = std::make_shared<arrow::StringBuilder>();
    fasta_h_builder = std::make_shared<arrow::StringBuilder>();
    rec_no_builder = std::make_shared<arrow::Int32Builder>();
    // rec_no_builder->Append(data.rec_no);
    // fasta_h_builder->Append(data.header);
    // fasta_seq_builder->Append(data.seq);
    // fasta_seq_builder->Finish(&seqArr);
    // fasta_h_builder->Finish(&hArr);
    // rec_no_builder->Finish(&recnoArr);
}

// Template
template <typename T>
std::shared_ptr<arrow::RecordBatch> ArrowWrapper::RecordBatchFrom(const T &str)
{
    std::shared_ptr<arrow::Array> seqArr, hArr, recnoArr; //}

    std::shared_ptr<arrow::Table> tab();
    // tab->AppendColumn("index", recnoArr);

    // std::shared_ptr<arrow::Buffer> recNoBuffer = arrow::Buffer::FromString(std::to_string(data.rec_no));
    // std::shared_ptr<arrow::Buffer> headerBuffer = arrow::Buffer::FromString(data.header);
    // std::shared_ptr<arrow::Buffer> seqBuffer = arrow::Buffer::FromString(data.seq);

    // std::shared_ptr<arrow::RecordBatch> rbatch = arrow::RecordBatch::Make(fasta_schema, 1, arrow::ArrayDataVector{entryArr});
    // omp_set_lock(&recBth_vecLock);
    // rbatch_vec->insert(rbatch_vec->end(), rbatch);
    // recBth_vec->insert(recBth_vec->end(), arrow::RecordBatch::Make(fasta_schema, 1, {recnoArr, hArr, seqArr}));
    // omp_unset_lock(&recBth_vecLock);
}

std::shared_ptr<arrow::Schema> ArrowWrapper::GetFASTASchema(void)
{
    return ArrowWrapper::GetSelf()->fasta_schema;
}

std::shared_ptr<arrow::KeyValueMetadata> ArrowWrapper::GetFASTAMetadata(void)
{
    return ArrowWrapper::GetSelf()->fasta_metadata;
}

std::vector<std::shared_ptr<arrow::Field>> ArrowWrapper::GetFASTAFields(void)
{
    std::vector<std::shared_ptr<arrow::Field>> fieldVector = {
        fasta_field0, fasta_field1, fasta_field2};
    return fieldVector;
}

arrow::ipc::IpcWriteOptions ArrowWrapper::GetArrowIPCOptions(void)
{
    return ArrowWrapper::GetSelf()->ipc_options;
}

int ArrowWrapper::GetColumnCount(const std::string_view &filename, char delim = '\t')
{
    std::ifstream file(filename.data());
    if (!file.is_open())
    {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return -1;
    }

    std::string line;
    if (std::getline(file, line))
    {
        std::stringstream ss(line);
        std::string column;
        int count = 0;
        while (std::getline(ss, column, delim))
        {
            count++;
        }
        return count;
    }
    else
    {
        std::cerr << "File is empty: " << filename << std::endl;
        return -1;
    }
}

std::shared_ptr<ArrowWrapper> ArrowWrapper::self = nullptr;