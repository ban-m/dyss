import cffi
import weakref
import os
import time
ffi = cffi.FFI()
ffi.cdef("""
typedef void* Dyss;
Dyss construct_dyss(const int num_scouts,const int num_packs, const char* ref_path, const char* model_pah,const char* param_path, int power,size_t querysize,size_t refsize);
int dyss_classify(Dyss d,int* query,size_t length);
int batch_classify(Dyss d,int** data, size_t* lengths,size_t num_of_data, int* result);
void dyss_destructor(Dyss d);
int is_null(Dyss d);
""")
if os.path.exists("./target/release/libdyss.so"):
    libdyss = ffi.dlopen("./target/release/libdyss.so")
elif os.path.exists("./target/release/libdyss.dylib"):
    libdyss = ffi.dlopen("./target/release/libdyss.dylib")
else:
    print("There are not dyss library. Please compile Dyss library first.")
    print("More specifically, after installing rust language and Cargo,:")
    print("cargo build --release")
    print("would be likely to solve the problem")
    sys.exit("Error:Couldn't open dyss library")

class Dyss(object):
    def __init__(self,
                 logger_name='Dyss',
                 num_scouts=14,
                 num_packs=3,
                 reference=b'./data/reference.fa',
                 model=b'./data/template_r9.4.model',
                 param=b'./data/parameters.csv',
                 power=9,
                 referencesize=200000):
        """
        A class to classify reads based on its raw signal.
        This class has two classifier method: dyss_classify and 
        batch_classify.
        Note that this class calls Rust function in the latter method,
        and the method call spawms as many process as possible.
        So it is not an appropriate way to create many thread and 
        run classify function on each thread.
        If you want to multithreading, please use the former method.
        The rust function itself static, i.e, it virtualy has no state.
        However, the thread safety is not garanteed.
        Wrapper class for Lock is strongly recommended.
        
        :param logger_name: the name to be used as a logger.
        :param num_scouts: the number of scuouts in extending procedure.
        :param num_packs: the number of packs in seeding procedure.
        :param reference: the reference to be used. First 'referencesize/2' base will be used.
        :param model: the model to be used to convert string into signal. HMM is used currently.
        :param power: the chunking power to be used in re-chunking procedure.
        :reference size: the total size of reference, i.e, template + reverse strand.
        """
        ref_path = ffi.new("char[]", reference.encode())
        model_path = ffi.new("char[]", model.encode())
        param_path = ffi.new("char[]", param.encode())
        self.count = 0
        self.positive = 0
        self.chunked = 0
        querysize = 250 # hard code
        self.classifier = libdyss.construct_dyss(num_scouts,
                                                 num_packs,
                                                 ref_path,
                                                 model_path,
                                                 param_path,
                                                 power,
                                                 querysize,
                                                 referencesize)
        # if libdyss.is_null(self.classifier) == 0:
        #     print("Dyss is null: Dyss doen't work.")
        #     print("Stop running script and review the settings immediately.")
        #     sys.exit("Error:Dyss is null")

    def batch_classify(self,queries):
        """ classify given queries parallely by using Rust language function.
        Note that the order of result is the same as input.
        Please do not shuffle the order in chunk_reads.
        param: queries: a python list of (channel_number,read_number,numpy int ndarray).
        """
        weakdict = weakref.WeakValueDictionary()
        length = len(queries)
        sent_data = ffi.new("int*[]",length)
        lengths = [0 for i in range(length)]
        for (i, (id,ch,r,x)) in enumerate(queries):
            lengths[i] = len(x)
            temp = ffi.new("int []",x)
            sent_data[i] = temp
            weakdict[temp] = sent_data
        sent_lengths = ffi.new("size_t[]",lengths)
        result = ffi.new("int[]",length)
        res = libdyss.batch_classify(self.classifier,sent_data,sent_lengths,length,result)
        self.count += length
        if res == 1: # worked correctly
            for x in result:
                if x == 1:
                    self.positive += 1
                elif x == 2:
                    self.chunked += 1
            return [(result[i],ch,r,id) for (i,(id,ch,r,data)) in enumerate(queries)]
        else:
            return None

    def free(self):
        libdyss.dyss_destructor(self.classifier)
        print("free dyss.Result:(positive,negative,chunked,total) = {}/{}/{}/{})".format(self.positive,self.count - self.chunked - self.positive,
                                                                                         self.chunked,self.count))

