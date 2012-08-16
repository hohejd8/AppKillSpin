#ifndef NO_FACTORY_LIMIT
bool ShouldInhibitFactoryRegistration(void (*const)()) {
  return true;
}
#endif
