
class BLASTWrapper
{
public:
    int num_threads = 4;
    std::string_view run_name;
    std::vector<std::string_view> queryFiles;
    std::vector<std::string_view> subjectFiles;
    std::string_view blastOptions = "-strand plus -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore score gaps frames qcovhsp sstrand qlen slen qseq sseq nident positive\"";

    // arrow::ArrayVector fArr{fHArr, fSArr};
    // std::shared_ptr<arrow::ChunkedArray> fArr_chunks = std::make_shared<arrow::ChunkedArray>(fArr);

private:
    int rec_count = 0;
    omp_lock_t rec_countLock;
    static std::shared_ptr<BLASTWrapper> self;
    BLASTWrapper();

public:
    ~BLASTWrapper();
    // void run_BLAST();
    // void (*PrintFastaBlock)(const FastaSequenceData &data, std::ostringstream *outputStream);

    void PrintFastaBlock(std::shared_ptr<FastaSequenceData> data, std::shared_ptr<std::ostringstream> outputStream);
    static void SetSelf(BLASTWrapper *_obj);
    static std::shared_ptr<BLASTWrapper> GetSelf(void);
    static BLASTWrapper *GetSelf_as_ObjPtr(void);
    std::shared_ptr<arrow::Schema> GetArrowFASTASchema(void);
    std::vector<std::shared_ptr<arrow::Field>> GetArrowFASTAFields(void);
    std::shared_ptr<arrow::KeyValueMetadata> GetArrowFASTAMetadata(void);
    std::shared_ptr<arrow::ipc::RecordBatchWriter> GetArrowRecordBatchWriter(std::shared_ptr<arrow::io::OutputStream> outputStream);
    arrow::ipc::IpcWriteOptions GetArrowIPCOptions(void);
    void SplitFASTA(std::string_view &fastaFile, int num_threads);

    // StreamFile() functions - Multithreaded file streaming, but NOT ordered
    // std::shared_ptr<T> StreamFile(const std::string_view &filename, const int num_threads, std::function<void(std::shared_ptr<FastaSequenceData>)> Entry_callback, const char *delim = "\n");
    template <typename T>
    std::shared_ptr<std::deque<std::shared_ptr<T>>> StreamFile(const std::string_view &filename, const int num_threads, std::function<void(std::shared_ptr<T>)> Entry_callback, const char *delim = "\n");
    std::shared_ptr<std::vector<std::shared_ptr<std::ostringstream>>> StreamFile(const std::string_view &filename, const int num_threads, const int &num_streams, const char *delim = "\n");
    std::shared_ptr<std::deque<std::shared_ptr<FastaSequenceData>>> StreamFile(const std::string_view &filename, const int num_threads, std::function<void(std::shared_ptr<FastaSequenceData>)> Entry_callback, const char *delim = ">");
    std::shared_ptr<arrow::RecordBatchVector> BLASTWrapper::StreamFile(const std::string_view &filename, const int num_threads, std::function<void(std::shared_ptr<arrow::RecordBatchVector>)> Entry_callback, const char *col_delim = "\t", const char *row_delim = "\n", const std::vector<std::string_view> &col_names = std::vector<std::string_view>());
    std::shared_ptr<BLASTHitData> BLASTWrapper::StreamFile(const std::string_view &filename, const int num_threads, std::function<void(std::shared_ptr<BLASTHitData>)> Entry_callback, const char *col_delim = "\t", const char *row_delim = "\n", const std::vector<std::string_view> &col_names = std::vector<std::string_view>());

private:
    // virtual void Init(void);
    // virtual int  Run(void);
    // virtual void Exit(void);
    std::shared_ptr<FastaSequenceData> stl_StreamFastaFile(const std::string_view &filename, const int num_threads, std::function<void(std::shared_ptr<FastaSequenceData>)> FASTAEntry_callback);
    void CreateBlastDatabase(const std::vector<std::string_view> &subjectFastaFiles);
    void PerformBlastQuery(const std::string_view &querySequence, const std::string_view &blastOptions);
    void LoadBlastDB(void);
};

BLASTWrapper::BLASTWrapper()
{

    SetSelf(this);
}

BLASTWrapper::~BLASTWrapper()
{
    std::cout << "~BLASTWrapper " << std::endl;

    // recBth_vec->clear();
    //  delete outputStream;
}

std::shared_ptr<BLASTWrapper> BLASTWrapper::GetSelf()
{
    if (BLASTWrapper::self.get() == nullptr)
    {
        BLASTWrapper *_obj = new BLASTWrapper();
        BLASTWrapper::SetSelf(_obj);
        // BLASTWrapper::self = static_cast<std::shared_ptr<BLASTWrapper>>(_obj);
    }
    return (BLASTWrapper::self);
}

BLASTWrapper *BLASTWrapper::GetSelf_as_ObjPtr(void)
{
    return BLASTWrapper::GetSelf().get();
}

void BLASTWrapper::SplitFASTA(std::string_view &fastaFile, int num_threads)
{
}

void BLASTWrapper::SetSelf(BLASTWrapper *_obj)
{
    BLASTWrapper::self = std::shared_ptr<BLASTWrapper>(_obj);
}

void BLASTWrapper::CreateBlastDatabase(const std::vector<std::string_view> &subjectFastaFiles)
{
    // ncbi::CSeqDB db(subjectFastaFiles, ncbi::CSeqDB::eNucleotide);
    // db.CreateNewBlastDb(subjectFastaFile.c_str(), eBlastEncodingNucleotide);
}

// Function to perform BLAST query for a single sequence
void BLASTWrapper::PerformBlastQuery(const std::string_view &querySequence, const std::string_view &blastOptions)
{
    // ncbi::blast::CBlast blast;
    //  Set BLAST options
    // blast.SetOptions(blastOptions.c_str());

    // Perform the query
    // ncbi::blast::CSearchResultSet results();
    // blast.Run(querySequence.c_str(), results);

    // Process the results as needed
    // ...
}

void BLASTWrapper::LoadBlastDB()
{
}

// Function to process a single FASTA block
void BLASTWrapper::PrintFastaBlock(std::shared_ptr<FastaSequenceData> data, std::shared_ptr<std::ostringstream> outputStream)
{
    if (outputStream != nullptr)
    {
        omp_set_lock(&outputStreamLock);
        // Perform processing on the FastaSequenceData object
        (*outputStream) << "No: " << data->rec_no << std::endl;
        (*outputStream) << "Header: " << data->header << std::endl;
        (*outputStream) << "Sequence: " << data->seq << std::endl;
        (*outputStream) << std::endl;
        outputStream->flush();
        omp_unset_lock(&outputStreamLock);
    }
}

std::shared_ptr<BLASTWrapper> BLASTWrapper::self = nullptr;