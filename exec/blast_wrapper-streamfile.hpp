template <typename T>
std::shared_ptr<std::deque<std::shared_ptr<T>>> BLASTWrapper::StreamFile(const std::string_view &filename, const int num_threads, std::function<void(std::shared_ptr<T>)> Entry_callback, const char *delim)
{
    return nullptr;
}

std::shared_ptr<std::vector<std::shared_ptr<std::ostringstream>>> BLASTWrapper::StreamFile(const std::string_view &filename, const int num_threads, const int &num_streams, const char *delim)
{
    return ArrowWrapper::GetSelf()->SplitFileIntoStreams(filename, num_threads, num_streams, delim);
}
// Function to stream FASTA file using mmap(), FASTAEntry_callback is a callback function to process each FASTA entry

std::shared_ptr<std::deque<std::shared_ptr<FastaSequenceData>>> BLASTWrapper::StreamFile(const std::string_view &filename, const int num_streams, std::function<void(std::shared_ptr<FastaSequenceData>)> Entry_callback, const char *delim)
{
    std::shared_ptr<std::deque<std::shared_ptr<FastaSequenceData>>> fastaPtr_q;
    omp_lock_t fastaPtr_qLock;

    fastaPtr_q = std::make_shared<std::deque<std::shared_ptr<FastaSequenceData>>>();
    fastaPtr_q->clear();

    // std::regex pattern(">([^\\n]+)\\n([ACTGU]+)"); //doesnt support multi-line
    std::regex pattern(">([^\\n]+)\\n((?:[ACTGU\\n]+)+)");

    omp_init_lock(&rec_countLock);
    omp_init_lock(&fastaPtr_qLock);

    std::shared_ptr<std::vector<std::shared_ptr<std::ostringstream>>> strStreams = StreamFile(filename, num_threads, num_streams, delim);
    omp_lock_t strStreamsLock;
    omp_init_lock(&strStreamsLock);
#pragma omp parallel shared(strStreams) shared(rec_count)
    {
#pragma omp for schedule(auto) nowait
        for (size_t i = 0; i < strStreams->size(); i++)
        {
            omp_set_lock(&strStreamsLock);
            std::shared_ptr<std::ostringstream> outputStream = strStreams->at(i);
            omp_unset_lock(&strStreamsLock);

            std::string_view fasta_entry = outputStream->str();
            if (!fasta_entry.empty())
            {
                std::smatch match;
                std::string tmp_str(fasta_entry);
                if (std::regex_search(tmp_str, match, pattern))
                {
                    omp_set_lock(&rec_countLock);
                    rec_count++;
                    omp_unset_lock(&rec_countLock);
                    FastaSequenceData data;
                    data.rec_no = rec_count;
                    data.header = match[1].str();
                    data.seq = match[2].str();
                    if (!data.header.empty() && !data.seq.empty())
                    {
                        Entry_callback(std::make_shared<FastaSequenceData>(data));
                    }
                }
            }
        }
#pragma omp barrier
    }

    std::cout << rec_count << std::endl;
    omp_destroy_lock(&rec_countLock);
    omp_destroy_lock(&fastaPtr_qLock);
    omp_destroy_lock(&strStreamsLock);

    return fastaPtr_q;
}

std::shared_ptr<arrow::RecordBatchVector> BLASTWrapper::StreamFile(const std::string_view &filename, const int num_threads, std::function<void(std::shared_ptr<arrow::RecordBatchVector>)> Entry_callback, const char *col_delim, const char *row_delim, const std::vector<std::string_view> &col_names)
{
    if (!col_names.empty())
    {
        // give a warning if number of columns is not equal to number of column names
    }

    std::shared_ptr<arrow::RecordBatchVector> rbatch_vec;
    omp_lock_t rbatch_vecLock;

    rbatch_vec = std::make_shared<arrow::RecordBatchVector>();
    rbatch_vec->clear();

    std::shared_ptr<std::vector<std::shared_ptr<std::ostringstream>>> strStreams = StreamFile(filename, num_threads, num_threads, row_delim);
    omp_lock_t strStreamsLock;
    omp_init_lock(&strStreamsLock);

    omp_init_lock(&rec_countLock);
    omp_init_lock(&rbatch_vecLock);
#pragma omp parallel shared(strStreams) shared(rec_count)
    {
#pragma omp for schedule(auto) nowait
        for (size_t i = 0; i < strStreams->size(); i++)
        {
            // split the entries into batches
            omp_set_lock(&strStreamsLock);
            std::shared_ptr<std::ostringstream> outputStream = strStreams->at(i);
            omp_unset_lock(&strStreamsLock);
            char *row_entries = strsep(*outputStream->str().c_str(), col_delim);
            std::cout << row_entries << std::endl;
        }

        omp_destroy_lock(&rec_countLock);
        omp_destroy_lock(&rbatch_vecLock);
        omp_destroy_lock(&strStreamsLock);
    }
#pragma omp barrier
}

std::shared_ptr<BLASTHitData> BLASTWrapper::StreamFile(const std::string_view &filename, const int num_threads, std::function<void(std::shared_ptr<BLASTHitData>)> Entry_callback, const char *col_delim, const char *row_delim, const std::vector<std::string_view> &col_names)
{
    if (!col_names.empty())
    {
        // give a warning if number of columns is not equal to number of column names
    }

    std::shared_ptr<arrow::RecordBatchVector> rbatch_vec;
    omp_lock_t rbatch_vecLock;

    rbatch_vec = std::make_shared<arrow::RecordBatchVector>();
    rbatch_vec->clear();

    std::shared_ptr<std::vector<std::shared_ptr<std::ostringstream>>> strStreams = StreamFile(filename, num_threads, num_threads, row_delim);
    omp_lock_t strStreamsLock;
    omp_init_lock(&strStreamsLock);

    omp_init_lock(&rec_countLock);
    omp_init_lock(&rbatch_vecLock);
#pragma omp parallel shared(strStreams) shared(rec_count)
    {
#pragma omp for schedule(auto) nowait
        for (size_t i = 0; i < strStreams->size(); i++)
        {
            // split the entries into batches
            omp_set_lock(&strStreamsLock);
            std::shared_ptr<std::ostringstream> outputStream = strStreams->at(i);
            omp_unset_lock(&strStreamsLock);
            char *row_entries = strsep(*outputStream->str().c_str(), col_delim);
            std::cout << row_entries << std::endl;
        }
#pragma omp barrier
    }
    omp_destroy_lock(&rec_countLock);
    omp_destroy_lock(&rbatch_vecLock);
    omp_destroy_lock(&strStreamsLock);
}