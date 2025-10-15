namespace ArcFrame.Solvers
{
    /// <summary>
    /// Reference data for implemented solutions. I don't wanna steal credit from anyone.
    /// </summary>
    public interface ISolverInfo
    {
        /// <summary>
        /// Reference key
        /// </summary>
        public string Key { get; }
        /// <summary>
        /// Main author
        /// </summary>
        public string Author { get; }
        /// <summary>
        /// DOI/URL etc
        /// </summary>
        public string Reference { get; }
    }
}
