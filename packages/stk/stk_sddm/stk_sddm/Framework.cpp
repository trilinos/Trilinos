//Add a dummy symbol so that linkers don't warn about this file being empty:
int dummy_framework_sddm_symbol(){return 0;}

// Dictionary &
//   static Dictionary &instance();

// Key *createKey(const KeyId &id, const std::string &name);
// Key *createTerminalKey(const KeyId &id, const std::string &name);

// Key *makeRoot(Key *key);


// Dictionary::instance()
// {
//   static Dictionary s_dictionary;

//   return s_dictionary;
// }


// Key *
// createKey(
//   const KeyId &         id,
//   const std::string &   name) 
// {
//   return Dictionary::instance().createKey(id, name);
// }


// Key *
// createTerminalKey(
//   const KeyId &         id,
//   const std::string &   name) 
// {
//   return Dictionary::instance().createTerminalKey(id, name);
// }


// Key *
// makeRoot(
//   Key *                 key) 
// {
//   return Dictionary::instance().makeRoot(key);
// }
