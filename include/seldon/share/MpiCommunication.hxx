#ifndef SELDON_FILE_MPI_COMMUNICATION_HXX

namespace Seldon
{
  
  const MPI::Datatype& GetMpiDataType(const Vector<bool>&);
  const MPI::Datatype& GetMpiDataType(const Vector<int>&);
  const MPI::Datatype& GetMpiDataType(const Vector<float>&);
  const MPI::Datatype& GetMpiDataType(const Vector<complex<float> >&);
  const MPI::Datatype& GetMpiDataType(const Vector<double>&);
  const MPI::Datatype& GetMpiDataType(const Vector<complex<double> >&);

  const MPI::Datatype& GetMpiDataType(const bool&);
  const MPI::Datatype& GetMpiDataType(const int&);
  const MPI::Datatype& GetMpiDataType(const float&);
  const MPI::Datatype& GetMpiDataType(const complex<float>&);
  const MPI::Datatype& GetMpiDataType(const double&);
  const MPI::Datatype& GetMpiDataType(const complex<double>&);
  
  template<class T>
  int GetRatioMpiDataType(const T&);
 
  template<class T>
  int GetRatioMpiDataType(const Vector<T>&);
 
  template<class T>
  int GetRatioMpiDataType(const complex<T>&);

  template<class T>
  int GetRatioMpiDataType(const Vector<complex<T> >&);

  template<class T>
  MPI::Request MpiIsend(const MPI::Comm& comm, T* x, Vector<int64_t>& xtmp,
                        int n, int proc, int tag);
  
  template<class T>
  MPI::Request MpiIsend(const MPI::Comm& comm, Vector<T>& x,
			Vector<int64_t>& xtmp,
                        int n, int proc, int tag);

  template<class T>
  MPI::Request MpiIrecv(const MPI::Comm& comm, T* x, Vector<int64_t>& xtmp,
                        int n, int proc, int tag);

  template<class T>
  MPI::Request MpiIrecv(const MPI::Comm& comm, Vector<T>& x,
			Vector<int64_t>& xtmp,
                        int n, int proc, int tag);
  
  template<class T>
  void MpiCompleteIrecv(T* x, Vector<int64_t>& xtmp, int n);
  
  template<class T>
  void MpiCompleteIrecv(Vector<T>& x, Vector<int64_t>& xtmp, int n);

  template<class T>
  void MpiSsend(const MPI::Comm& comm, T* x, Vector<int64_t>& xtmp,
                int n, int proc, int tag);
  
  template<class T>
  void MpiSsend(const MPI::Comm& comm, Vector<T>& x, Vector<int64_t>& xtmp,
                int n, int proc, int tag);
                       
  template<class T>
  void MpiSend(const MPI::Comm& comm, T* x, Vector<int64_t>& xtmp,
               int n, int proc, int tag);
  
  template<class T>
  void MpiSend(const MPI::Comm& comm, Vector<T>& x, Vector<int64_t>& xtmp,
               int n, int proc, int tag);
  
  template<class T>
  void MpiGather(const MPI::Comm& comm, T* x, Vector<int64_t>& xtmp,
                 T* y, int n, int proc);
  
  template<class T>
  void MpiGather(const MPI::Comm& comm, Vector<T>& x, Vector<int64_t>& xtmp,
                 Vector<T>& y, int n, int proc);
  
  template<class T>
  void MpiAllreduce(const MPI::Comm& comm, T* x, Vector<int64_t>& xtmp,
                    T* y, int n, const MPI::Op& op);
  
  template<class T>
  void MpiAllreduce(const MPI::Comm& comm, Vector<T>& x,
		    Vector<int64_t>& xtmp,
                    Vector<T>& y, int n, const MPI::Op& op);
  
  template<class T>
  void MpiReduce(const MPI::Comm& comm, T* x, Vector<int64_t>& xtmp,
                 T* y, int n, const MPI::Op& op, int proc);
  
  template<class T>
  void MpiReduce(const MPI::Comm& comm, Vector<T>& x, Vector<int64_t>& xtmp,
                 Vector<T>& y, int n, const MPI::Op& op, int proc);
  
  template<class T>
  void MpiRecv(const MPI::Comm& comm, T* x, Vector<int64_t>& xtmp,
               int n, int proc, int tag, MPI::Status& status);
  
  template<class T>
  void MpiRecv(const MPI::Comm& comm, Vector<T>& x, Vector<int64_t>& xtmp,
               int n, int proc, int tag, MPI::Status& status);

  template<class T>
  void MpiBcast(const MPI::Comm& comm, T* x,
		Vector<int64_t>& xtmp, int n, int proc);
  
  template<class T>
  void MpiBcast(const MPI::Comm& comm, Vector<T>& x,
		Vector<int64_t>& xtmp, int n, int proc);
  
}

#define SELDON_FILE_MPI_COMMUNICATION_HXX
#endif
