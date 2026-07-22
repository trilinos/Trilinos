// Copyright(C) 2025 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <fmt/core.h>
#include <fmt/ostream.h>
#include <getopt.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cctype>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include <aws/core/Aws.h>
#include <aws/core/auth/AWSCredentialsProvider.h>
#include <aws/core/auth/AWSCredentialsProviderChain.h>
#include <aws/core/utils/logging/LogLevel.h>
#include <aws/core/utils/stream/PreallocatedStreamBuf.h>

#include <aws/s3/S3Client.h>

#include <aws/s3/model/CreateBucketRequest.h>
#include <aws/s3/model/DeleteBucketPolicyRequest.h>
#include <aws/s3/model/DeleteBucketRequest.h>
#include <aws/s3/model/DeleteObjectRequest.h>
#include <aws/s3/model/GetBucketPolicyRequest.h>
#include <aws/s3/model/GetObjectRequest.h>
#include <aws/s3/model/HeadBucketRequest.h>
#include <aws/s3/model/ListObjectsRequest.h>
#include <aws/s3/model/PutBucketPolicyRequest.h>
#include <aws/s3/model/PutObjectRequest.h>

#include <aws/core/utils/threading/Executor.h>
#include <aws/transfer/TransferHandle.h>
#include <aws/transfer/TransferManager.h>

#include <Ioss_ParallelUtils.h>
#include <Ioss_Property.h>
#include <Ioss_PropertyManager.h>
#include <Ioss_Utils.h>

#include "Ios3_AwsHelpers.h"

namespace Ios3 {
  namespace helpers {

    const std::string env_base_name{"IOSS_S3_"};
    const std::string env_name_endpoint{env_base_name + "ENDPOINT"};
    const std::string env_name_profile{env_base_name + "PROFILE"};
    const std::string env_name_ca_file{env_base_name + "CA_FILE"};
    const std::string env_name_use_ca_file{env_base_name + "USE_CA_FILE"};
    const std::string env_name_use_transfer_manager{env_base_name + "USE_TRANSFER_MANAGER"};
    const std::string env_name_enable_aws_tracing{env_base_name + "ENABLE_AWS_TRACING"};
    const std::string env_name_disable_ec2_lookup{env_base_name + "DISABLE_EC2_LOOKUP"};

    const std::string prop_base_name{"S3_"};
    const std::string prop_name_endpoint{prop_base_name + "ENDPOINT"};
    const std::string prop_name_profile{prop_base_name + "PROFILE"};
    const std::string prop_name_ca_file{prop_base_name + "CA_FILE"};
    const std::string prop_name_use_ca_file{prop_base_name + "USE_CA_FILE"};
    const std::string prop_name_use_transfer_manager{prop_base_name + "USE_TRANSFER_MANAGER"};
    const std::string prop_name_enable_aws_tracing{prop_base_name + "ENABLE_AWS_TRACING"};
    const std::string prop_name_disable_ec2_lookup{prop_base_name + "DISABLE_EC2_LOOKUP"};

    void print_params(const HelperParameters &params)
    {
      fmt::print(stdout, "INFO: endpoint             = {}\n", params.endpoint);
      fmt::print(stdout, "INFO: profile              = {}\n", params.profile);
      fmt::print(stdout, "INFO: use_ca_file          = {}\n", params.use_ca_file);
      fmt::print(stdout, "INFO: ca_file              = {}\n", params.ca_file);
      fmt::print(stdout, "INFO: use_transfer_manager = {}\n", params.use_transfer_manager);
      fmt::print(stdout, "INFO: enable_aws_tracing   = {}\n", params.enable_aws_tracing);
      fmt::print(stdout, "INFO: disable_ec2_lookup   = {}\n", params.disable_ec2_lookup);
    }

    static const char AwsHelperCredentialsProviderChainTag[] = "AwsHelperCredentialsProviderChain";

    /**
     * Creates an AWSCredentialsProviderChain which only uses EnvironmentAWSCredentialsProvider
     * and ProfileConfigFileAWSCredentialsProvider in that order.
     */
    class AwsHelperCredentialsProviderChain : public Aws::Auth::AWSCredentialsProviderChain
    {
    public:
      /**
       * Initializes the provider chain with EnvironmentAWSCredentialsProvider
       * and ProfileConfigFileAWSCredentialsProvider in that order.
       */
      AwsHelperCredentialsProviderChain() : AWSCredentialsProviderChain()
      {
        AddProvider(Aws::MakeShared<Aws::Auth::EnvironmentAWSCredentialsProvider>(
            AwsHelperCredentialsProviderChainTag));
        AddProvider(Aws::MakeShared<Aws::Auth::ProfileConfigFileAWSCredentialsProvider>(
            AwsHelperCredentialsProviderChainTag));
      }
      /**
       * Initializes the provider chain with EnvironmentAWSCredentialsProvider
       * and ProfileConfigFileAWSCredentialsProvider in that order.
       */
      AwsHelperCredentialsProviderChain(const char *profile) : AWSCredentialsProviderChain()
      {
        AddProvider(Aws::MakeShared<Aws::Auth::EnvironmentAWSCredentialsProvider>(
            AwsHelperCredentialsProviderChainTag));
        AddProvider(Aws::MakeShared<Aws::Auth::ProfileConfigFileAWSCredentialsProvider>(
            AwsHelperCredentialsProviderChainTag, profile));
      }
    };

    static int             context_count = 0;
    static Aws::SDKOptions options;

    void getPropertiesFromEnvVars(Ioss::PropertyManager     &properties,
                                  const Ioss::ParallelUtils &utils)
    {
      std::string env_value;
      if (utils.get_environment(env_name_endpoint, env_value, utils.parallel_size() > 1)) {
        properties.add(Ioss::Property(prop_name_endpoint, env_value));
      }
      if (utils.get_environment(env_name_profile, env_value, utils.parallel_size() > 1)) {
        properties.add(Ioss::Property(prop_name_profile, env_value));
      }
      if (utils.get_environment(env_name_ca_file, env_value, utils.parallel_size() > 1)) {
        properties.add(Ioss::Property(prop_name_ca_file, env_value));
      }
      if (utils.get_environment(env_name_use_ca_file, env_value, utils.parallel_size() > 1)) {
        int         prop_value   = 0;
        std::string up_env_value = Ioss::Utils::uppercase(env_value);
        if (up_env_value == "1" || up_env_value == "TRUE" || up_env_value == "YES" ||
            up_env_value == "ON") {
          prop_value = 1;
        }
        properties.add(Ioss::Property(prop_name_use_ca_file, prop_value));
      }
      if (utils.get_environment(env_name_use_transfer_manager, env_value,
                                utils.parallel_size() > 1)) {
        int         prop_value   = 0;
        std::string up_env_value = Ioss::Utils::uppercase(env_value);
        if (up_env_value == "1" || up_env_value == "TRUE" || up_env_value == "YES" ||
            up_env_value == "ON") {
          prop_value = 1;
        }
        properties.add(Ioss::Property(prop_name_use_transfer_manager, prop_value));
      }
      if (utils.get_environment(env_name_enable_aws_tracing, env_value,
                                utils.parallel_size() > 1)) {
        int         prop_value   = 0;
        std::string up_env_value = Ioss::Utils::uppercase(env_value);
        if (up_env_value == "1" || up_env_value == "TRUE" || up_env_value == "YES" ||
            up_env_value == "ON") {
          prop_value = 1;
        }
        properties.add(Ioss::Property(prop_name_enable_aws_tracing, prop_value));
      }
      if (utils.get_environment(env_name_disable_ec2_lookup, env_value,
                                utils.parallel_size() > 1)) {
        int         prop_value   = 1;
        std::string up_env_value = Ioss::Utils::uppercase(env_value);
        if (up_env_value != "1" && up_env_value != "TRUE" && up_env_value != "YES" &&
            up_env_value != "ON") {
          prop_value = 0;
        }
        properties.add(Ioss::Property(prop_name_disable_ec2_lookup, prop_value));
      }
    }

    void getParamsFromProperties(Ioss::PropertyManager &properties, HelperParameters &params)
    {
      if (properties.exists(prop_name_endpoint)) {
        params.endpoint = properties.get(prop_name_endpoint).get_string();
      }
      if (properties.exists(prop_name_profile)) {
        params.profile = properties.get(prop_name_profile).get_string();
      }
      if (properties.exists(prop_name_ca_file)) {
        params.ca_file = properties.get(prop_name_ca_file).get_string();
      }
      if (properties.exists(prop_name_use_ca_file)) {
        params.use_ca_file = (1 == properties.get(prop_name_use_ca_file).get_int());
      }
      if (properties.exists(prop_name_use_transfer_manager)) {
        params.use_transfer_manager =
            (1 == properties.get(prop_name_use_transfer_manager).get_int());
      }
      if (properties.exists(prop_name_enable_aws_tracing)) {
        params.enable_aws_tracing = (1 == properties.get(prop_name_enable_aws_tracing).get_int());
      }
      if (properties.exists(prop_name_disable_ec2_lookup)) {
        params.disable_ec2_lookup = (1 == properties.get(prop_name_disable_ec2_lookup).get_int());
      }
    }

    std::shared_ptr<HelperContext> createContext(const HelperParameters &params)
    {
      std::shared_ptr<HelperContext> context = std::make_shared<HelperContext>();

      context->use_transfer_manager = params.use_transfer_manager;

      if (params.disable_ec2_lookup) {
        // AWS tries to call out to an EC2 server.  If your S3 service doesn't
        // have EC2, setting this envvar will eliminate a 1 second
        // startup delay while it waits for the connection to timeout.
        setenv("AWS_EC2_METADATA_DISABLED", "true", 1 /* overwrite */);
      }

      // this is not thread safe
      if (++context_count == 1) {
        // If the object store closes our connections, a sigpipe is
        // generated inside Curl that isn't handled by Curl.  This
        // option installs an AWS sigpipe handler that prevents an
        // exit due to the unhandled signal.
        options.httpOptions.installSigPipeHandler = true;
        if (params.enable_aws_tracing) {
          options.loggingOptions.logLevel = Aws::Utils::Logging::LogLevel::Trace;
        }
        InitAPI(options);
      }

      Aws::Client::ClientConfiguration config;
      config.endpointOverride = params.endpoint;
      config.profileName      = params.profile;
      config.requestTimeoutMs = 100000;
      if (params.use_ca_file) {
        config.caFile = params.ca_file;
      }

      // This creates a credentials provider chain that first looks for
      // the AWS_ACCESS_KEY_ID and AWS_SECRET_ACCESS_KEY environment
      // variables and then looks for a specific profile from the user's
      // config files ($HOME/.aws/credentials).  The default provider is a
      // provider chain that instantiates all the potential providers with
      // the default constructor.  It doesn't seem to use the
      // ClientConfiguration parameters, so it always looks for the
      // "default" profile.  This allows us to specify the provider we
      // want and any parameters we want.
      auto credProvider = Aws::MakeShared<AwsHelperCredentialsProviderChain>(
          "CredProvider", params.profile.c_str());

      context->client = Aws::MakeShared<Aws::S3::S3Client>(
          "S3Client", credProvider, config,
          Aws::Client::AWSAuthV4Signer::PayloadSigningPolicy::Never,
          false /* disable virtual addressing - otherwise client requests timeout */);

      context->executor =
          Aws::MakeShared<Aws::Utils::Threading::PooledThreadExecutor>("executor", 25);

      return context;
    }

    void destroyContext(std::shared_ptr<HelperContext> context)
    {
      context->transfer_manager.reset();
      context->client.reset();
      ;

      // this is not thread safe
      if (--context_count == 0) {
        ShutdownAPI(options);
      }

      // destroy the context held by the shared pointer
      context.reset();
    }

    namespace {

      // derived from example found here: https://en.cppreference.com/w/cpp/string/byte/tolower
      std::string tolower(std::string s)
      {
        std::transform(s.begin(), s.end(), s.begin(),
                       [](unsigned char c) { return std::tolower(c); });
        return s;
      }

      std::string replace_aws_illegal_chars(std::string s)
      {
        std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
          char legal = c;
          if (!std::isalpha(c) && !std::isdigit(c) && c != '-') {
            legal = '-';
          }
          return legal;
        });
        return s;
      }

    } // namespace

    std::string cleanBucketName(const std::string &name)
    {
      std::string trimmed_name{name};
      auto        pos = name.rfind("/");
      if (pos != std::string::npos) {
        trimmed_name = name.substr(pos + 1);
      }

      std::string::size_type n     = 0;
      std::string            valid = trimmed_name;

      valid = tolower(valid);

      // '.' is discouraged so replace with '-'
      n = valid.find('.');
      while (n != std::string::npos) {
        valid.at(n) = '-';
        n           = valid.find('.', n);
      }
      // can't start with "xn--"
      if (valid.find("xn--") == 0) {
        valid = 'x' + valid;
      }
      // can't end with "--s3alias"
      n = valid.rfind("--s3alias");
      if (n + 9 == valid.length()) {
        valid = valid.substr(0, n);
      }
      // can't end with "--ol-s3"
      n = valid.rfind("--ol-s3");
      if (n + 7 == valid.length()) {
        valid = valid.substr(0, n);
      }

      if (valid.length() < 3) {
        // minimum bucket name length is 3 chars
        valid = valid + "-zzz";
      }
      else if (valid.length() > 63) {
        // maximum bucket name length is 63 chars
        valid = valid.substr(0, 63);
      }

      // must start with a letter or number
      if (!std::isalnum(valid.front())) {
        if (valid.length() < 63) {
          valid = '0' + valid;
        }
        else {
          valid.front() = '0';
        }
      }
      // must end with a letter or number
      if (!std::isalnum(valid.back())) {
        if (valid.length() < 63) {
          valid = valid + '9';
        }
        else {
          valid.back() = '9';
        }
      }

      valid = replace_aws_illegal_chars(valid);

      return valid;
    }

    bool createBucket(std::shared_ptr<HelperContext> context, const std::string &bucket)
    {
      bool success = true;

      std::string clean_bucket_name = cleanBucketName(bucket);

      Aws::S3::Model::CreateBucketRequest request;
      request.SetBucket(clean_bucket_name);

      Aws::S3::Model::CreateBucketOutcome outcome = context->client->CreateBucket(request);
      if (!outcome.IsSuccess()) {
        auto err = outcome.GetError();
        fmt::print(stderr, "Error: CreateBucket: {}: {}", err.GetExceptionName(), err.GetMessage());
        success = false;
      }

      return success;
    }

    bool waitBucket(std::shared_ptr<HelperContext> context, const std::string &bucket,
                    uint64_t wait_usec)
    {
      bool success = false;

      std::string clean_bucket_name = cleanBucketName(bucket);

      // convert microseconrds to nanoseconds then divide into 10 iterations
      uint64_t per_sleep_nsec = (wait_usec * 1000) / 10;
      unsigned timeoutCount   = 0;
      while (timeoutCount++ < 10) {
        Aws::S3::Model::HeadBucketRequest headBucketRequest;
        headBucketRequest.SetBucket(clean_bucket_name);
        Aws::S3::Model::HeadBucketOutcome headBucketOutcome =
            context->client->HeadBucket(headBucketRequest);
        if (headBucketOutcome.IsSuccess()) {
          success = true;
          break;
        }
        std::this_thread::sleep_for(std::chrono::nanoseconds(per_sleep_nsec));
      }

      return success;
    }

    bool deleteBucket(std::shared_ptr<HelperContext> context, const std::string &bucket)
    {
      bool success = true;

      std::string clean_bucket_name = cleanBucketName(bucket);

      Aws::S3::Model::DeleteBucketRequest request;
      request.SetBucket(clean_bucket_name);

      Aws::S3::Model::DeleteBucketOutcome outcome = context->client->DeleteBucket(request);
      if (!outcome.IsSuccess()) {
        auto err = outcome.GetError();
        fmt::print(stderr, "Error: DeleteBucket: {}: {}", err.GetExceptionName(), err.GetMessage());
        success = false;
      }

      return success;
    }

    bool listBuckets(std::shared_ptr<HelperContext> context, std::vector<std::string> &bucket_names)
    {
      bool success = true;

      auto outcome = context->client->ListBuckets();
      if (outcome.IsSuccess()) {
        auto buckets = outcome.GetResult().GetBuckets();
        for (auto &&b : buckets) {
          bucket_names.push_back(b.GetName());
        }
      }
      else {
        auto err = outcome.GetError();
        fmt::print(stderr, "Error: ListBuckets: {}: {}", err.GetExceptionName(), err.GetMessage());
        success = false;
      }

      return success;
    }

    template <typename T, typename A>
    bool putValue(std::shared_ptr<HelperContext> context, const std::string &bucket,
                  const std::string &key, const std::vector<T, A> &value)
    {
      bool success = true;

      // PreallocatedStreamBuf doesn't have a read-only behavior, so the underlying data can't be
      // const.
      std::vector<T, A>                        *vec = const_cast<std::vector<T, A> *>(&value);
      Aws::Utils::Stream::PreallocatedStreamBuf streamBuffer(vec->data(), vec->size());
      std::shared_ptr<Aws::IOStream>            input_stream =
          Aws::MakeShared<Aws::IOStream>("SomeTag", &streamBuffer);

      std::string clean_bucket_name = cleanBucketName(bucket);

      if (context->use_transfer_manager) {
        Aws::Transfer::TransferManagerConfiguration transfer_config(context->executor.get());
        const uint64_t                              MB5 = 5 * 1024 * 1024UL;
        transfer_config.bufferSize                      = MB5;
        transfer_config.s3Client                        = context->client;
        context->transfer_manager = Aws::Transfer::TransferManager::Create(transfer_config);

        std::shared_ptr<Aws::Transfer::TransferHandle> uploadHandle =
            context->transfer_manager->UploadFile(input_stream, clean_bucket_name, key,
                                                  "binary/octet-stream",
                                                  Aws::Map<Aws::String, Aws::String>());
        uploadHandle->WaitUntilFinished();

        if (uploadHandle->GetStatus() != Aws::Transfer::TransferStatus::COMPLETED) {
          auto err = uploadHandle->GetLastError();
          fmt::print(stderr, "Error: TransferManager::UploadFile: {}: {}", err.GetExceptionName(),
                     err.GetMessage());
          success = false;
        }
      }
      else {
        Aws::S3::Model::PutObjectRequest request;
        request.SetBucket(clean_bucket_name);
        request.SetKey(key);

        request.SetBody(input_stream);
        request.SetContentLength(value.size());

        Aws::S3::Model::PutObjectOutcome outcome = context->client->PutObject(request);
        if (!outcome.IsSuccess()) {
          auto err = outcome.GetError();
          fmt::print(stderr, "Error: PutObject: {}: {}", err.GetExceptionName(), err.GetMessage());
          success = false;
        }
      }

      return success;
    }

    template <typename T, typename A>
    bool getValue(std::shared_ptr<HelperContext> context, const std::string &bucket,
                  const std::string &key, std::vector<T, A> &value)
    {
      bool success = true;

      std::string clean_bucket_name = cleanBucketName(bucket);

      Aws::S3::Model::GetObjectRequest request;
      request.SetBucket(clean_bucket_name);
      request.SetKey(key);

      Aws::S3::Model::GetObjectOutcome outcome = context->client->GetObject(request);
      if (outcome.IsSuccess()) {
        auto   &retrieved_obj_body = outcome.GetResultWithOwnership().GetBody();
        int64_t contentLength      = outcome.GetResultWithOwnership().GetContentLength();
        int64_t total_bytes_read   = 0;
        int64_t bytes_left         = contentLength;
        if (contentLength > (int64_t)value.capacity()) {
          value.reserve(contentLength);
          value.resize(contentLength);
        }
        while ((bytes_left > 0) && (total_bytes_read < (int64_t)value.capacity())) {
          int64_t bytes_read =
              retrieved_obj_body.readsome((char *)&value[total_bytes_read], bytes_left);
          total_bytes_read += bytes_read;
          bytes_left -= bytes_read;
        }
      }
      else {
        auto err = outcome.GetError();
        fmt::print(stderr, "Error: GetObject: {}: {}", err.GetExceptionName(), err.GetMessage());
        success = false;
      }

      return success;
    }

    bool putValue(std::shared_ptr<HelperContext> context, const std::string &bucket,
                  const std::string &key, const std::string &filename)
    {
      bool success = true;

      std::string clean_bucket_name = cleanBucketName(bucket);

      if (context->use_transfer_manager) {
        Aws::Transfer::TransferManagerConfiguration transfer_config(context->executor.get());
        const uint64_t                              GB5 = 5 * 1024 * 1024 * 1024UL;
        transfer_config.bufferSize                      = GB5;
        transfer_config.s3Client                        = context->client;
        context->transfer_manager = Aws::Transfer::TransferManager::Create(transfer_config);

        std::shared_ptr<Aws::Transfer::TransferHandle> uploadHandle =
            context->transfer_manager->UploadFile(filename, clean_bucket_name, key,
                                                  "binary/octet-stream",
                                                  Aws::Map<Aws::String, Aws::String>());
        uploadHandle->WaitUntilFinished();

        if (uploadHandle->GetStatus() != Aws::Transfer::TransferStatus::COMPLETED) {
          auto err = uploadHandle->GetLastError();
          fmt::print(stderr, "Error: TransferManager::UploadFile: {}: {}", err.GetExceptionName(),
                     err.GetMessage());
          success = false;
        }
      }
      else {
        fmt::print(stderr,
                   "Error: PutValue: Putting from a file only works with the transfer manager.");
        success = false;
      }

      return success;
    }

    bool getValue(std::shared_ptr<HelperContext> context, const std::string &bucket,
                  const std::string &key, std::string &filename)
    {
      bool success = true;

      std::string clean_bucket_name = cleanBucketName(bucket);

      Aws::S3::Model::GetObjectRequest request;
      request.SetBucket(clean_bucket_name);
      request.SetKey(key);

      Aws::S3::Model::GetObjectOutcome outcome = context->client->GetObject(request);
      if (outcome.IsSuccess()) {
        std::ofstream   obj_ofs(filename);
        auto           &retrieved_obj_body    = outcome.GetResultWithOwnership().GetBody();
        constexpr int   max_read              = 16385;
        char            object_data[max_read] = {0};
        int64_t         contentLength         = outcome.GetResultWithOwnership().GetContentLength();
        std::streamsize bytes_left            = contentLength;

        while (bytes_left > 0) {
          std::streamsize bytes_read = retrieved_obj_body.readsome(object_data, max_read);
          bytes_left -= bytes_read;
          obj_ofs.write(object_data, bytes_read);
        }
      }
      else {
        auto err = outcome.GetError();
        fmt::print(stderr, "Error: GetObject: {}: {}", err.GetExceptionName(), err.GetMessage());
        success = false;
      }

      return success;
    }

    bool deleteValue(std::shared_ptr<HelperContext> context, const std::string &bucket,
                     const std::string &key)
    {
      bool success = true;

      std::string clean_bucket_name = cleanBucketName(bucket);

      Aws::S3::Model::DeleteObjectRequest request;
      request.SetBucket(clean_bucket_name);
      request.SetKey(key);

      Aws::S3::Model::DeleteObjectOutcome outcome = context->client->DeleteObject(request);
      if (!outcome.IsSuccess()) {
        auto err = outcome.GetError();
        fmt::print(stderr, "Error: DeleteObject: {}: {}", err.GetExceptionName(), err.GetMessage());
        success = false;
      }

      return success;
    }

    bool listKeys(std::shared_ptr<HelperContext> context, const std::string &bucket,
                  std::vector<std::string> &keys)
    {
      bool success = true;

      std::string clean_bucket_name = cleanBucketName(bucket);

      Aws::S3::Model::ListObjectsRequest request;
      request.SetBucket(clean_bucket_name);

      Aws::S3::Model::ListObjectsOutcome outcome = context->client->ListObjects(request);
      if (outcome.IsSuccess()) {
        for (auto &&b : outcome.GetResult().GetContents()) {
          keys.push_back(b.GetKey());
        }
      }
      else {
        auto err = outcome.GetError();
        fmt::print(stderr, "Error: ListObjects: {}: {}", err.GetExceptionName(), err.GetMessage());
        success = false;
      }

      return success;
    }

    bool listKeys(std::shared_ptr<HelperContext> context, const std::string &bucket,
                  const std::string &key_prefix, std::vector<std::string> &keys)
    {
      bool success = true;

      std::string clean_bucket_name = cleanBucketName(bucket);

      Aws::S3::Model::ListObjectsRequest request;
      request.SetBucket(clean_bucket_name);
      request.SetPrefix(key_prefix);

      Aws::S3::Model::ListObjectsOutcome outcome = context->client->ListObjects(request);
      if (outcome.IsSuccess()) {
        for (auto &&b : outcome.GetResult().GetContents()) {
          keys.push_back(b.GetKey());
        }
      }
      else {
        auto err = outcome.GetError();
        fmt::print(stderr, "Error: ListObjects: {}: {}", err.GetExceptionName(), err.GetMessage());
        success = false;
      }

      return success;
    }

    bool putBucketPolicy(std::shared_ptr<HelperContext> context, const std::string &bucket,
                         const std::string &policy)
    {
      bool success = true;

      std::string clean_bucket_name = cleanBucketName(bucket);

      std::shared_ptr<Aws::StringStream> request_body = Aws::MakeShared<Aws::StringStream>("");
      *request_body << policy;

      Aws::S3::Model::PutBucketPolicyRequest request;
      request.SetBucket(clean_bucket_name);
      request.SetBody(request_body);

      Aws::S3::Model::PutBucketPolicyOutcome outcome = context->client->PutBucketPolicy(request);
      if (!outcome.IsSuccess()) {
        auto err = outcome.GetError();
        fmt::print(stderr, "Error: PutBucketPolicy: {}: {}", err.GetExceptionName(),
                   err.GetMessage());
        success = false;
      }

      return success;
    }

    bool getBucketPolicy(std::shared_ptr<HelperContext> context, const std::string &bucket,
                         std::string &policy)
    {
      bool success = true;

      std::string clean_bucket_name = cleanBucketName(bucket);

      Aws::S3::Model::GetBucketPolicyRequest request;
      request.SetBucket(clean_bucket_name);

      Aws::S3::Model::GetBucketPolicyOutcome outcome = context->client->GetBucketPolicy(request);
      if (outcome.IsSuccess()) {
        outcome.GetResult().GetPolicy() >> policy;
      }
      else {
        auto err = outcome.GetError();
        fmt::print(stderr, "Error: GetBucketPolicy: {}: {}", err.GetExceptionName(),
                   err.GetMessage());
        success = false;
      }

      return success;
    }

    bool deleteBucketPolicy(std::shared_ptr<HelperContext> context, const std::string &bucket)
    {
      bool success = true;

      std::string clean_bucket_name = cleanBucketName(bucket);

      Aws::S3::Model::DeleteBucketPolicyRequest request;
      request.SetBucket(clean_bucket_name);

      Aws::S3::Model::DeleteBucketPolicyOutcome outcome =
          context->client->DeleteBucketPolicy(request);
      if (!outcome.IsSuccess()) {
        auto err = outcome.GetError();
        fmt::print(stderr, "Error: DeleteBucketPolicy: {}: {}", err.GetExceptionName(),
                   err.GetMessage());
        success = false;
      }

      return success;
    }

    /*
     * Explicit Template Instantiation
     */

    template bool putValue(std::shared_ptr<HelperContext> context, const std::string &bucket,
                           const std::string                                               &key,
                           const std::vector<unsigned char, std::allocator<unsigned char>> &value);
    template bool putValue(std::shared_ptr<HelperContext> context, const std::string &bucket,
                           const std::string &key, const UninitializedVector<unsigned char> &value);

    template bool getValue(std::shared_ptr<HelperContext> context, const std::string &bucket,
                           const std::string                                         &key,
                           std::vector<unsigned char, std::allocator<unsigned char>> &value);
    template bool getValue(std::shared_ptr<HelperContext> context, const std::string &bucket,
                           const std::string &key, UninitializedVector<unsigned char> &value);

  } // namespace helpers
} // namespace Ios3
